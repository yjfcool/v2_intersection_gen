#include "connectivity_generator.h"
#include "curve/hermite_init.h"
#include "curve/curve_utils.h"
#include "constraints/fence_check.h"
#include "constraints/intersection_check.h"
#include <chrono>
#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>

#include "iodata_shapefile.h"
#include "shapefile.hpp"

// ── IntersectionInput helpers ─────────────────────────────────
const Lane* IntersectionInput::findLane(const LaneId&id)const{for(auto&l:lanes)if(l.id==id)return&l;return nullptr;}
const LaneGroup* IntersectionInput::findGroup(const LaneGroupId&id)const{for(auto&g:lane_groups)if(g.id==id)return&g;return nullptr;}
const LaneEdge* IntersectionInput::findEdge(const LaneEdgeId&id)const{for(auto&e:lane_edges)if(e.id==id)return&e;return nullptr;}
bool IntersectionInput::laneGroupExists(const LaneGroupId&id)const{return findGroup(id)!=nullptr;}

static const LaneGroup* groupOf(const IntersectionInput&inp,const LaneId&lid,GroupRole role){
    for(auto&g:inp.lane_groups){if(g.role!=role)continue;for(auto&id:g.lanes)if(id==lid)return&g;}
    return nullptr;
}
// Derive endpoint & tangent direction from Lane::geometry
// Entry lane geometry: outside → intersection edge (last point = entry point into junction)
// Exit  lane geometry: intersection edge → outside (first point = exit point from junction)
std::pair<Vec2d,Vec2d> IntersectionInput::entryPtDir(const LaneId&lid)const{
    auto*l=findLane(lid);
    if(!l||l->geometry.points.empty()){
        std::cout << "[WARN] entrylane:" << lid << " no geometry !" << std::endl;
        // fallback to group ref_point if no geometry available
        auto*g=groupOf(*this,lid,GroupRole::Entry);
        return //g?std::make_pair(g->ref_point,g->direction) :
                std::make_pair(Vec2d(0,0),Vec2d(1,0));
    }
    return{laneEntryPoint(l->geometry), laneEntryTangent(l->geometry)};
}
std::pair<Vec2d,Vec2d> IntersectionInput::exitPtDir(const LaneId&lid)const{
    auto*l=findLane(lid);
    if(!l||l->geometry.points.empty()){
        std::cout << "[WARN] exitlane:" << lid << " no geometry !" << std::endl;
        auto*g=groupOf(*this,lid,GroupRole::Exit);
        return //g?std::make_pair(g->ref_point,g->direction) :
                std::make_pair(Vec2d(10,0),Vec2d(1,0));
    }
    return{laneExitPoint(l->geometry), laneExitTangent(l->geometry)};
}

// ── GlobalCoordinator ─────────────────────────────────────────
int GlobalCoordinator::turnPriority(TurnType t){
    switch(t){case TurnType::Straight:return 0;case TurnType::TurnLeft:return 1;
    case TurnType::TurnRight:return 2;default:return 3;}
}
void GlobalCoordinator::build(const std::vector<Connectivity>&conns,const IntersectionInput&){
    std::map<int,OptGroup>pm;
    for(auto&c:conns){int p=turnPriority(c.turn_type);pm[p].conn_ids.push_back(c.id);pm[p].priority=p;}
    groups_.clear();
    for(auto&[p,g]:pm)groups_.push_back(g);
}
void GlobalCoordinator::addSoftObstacles(SDFField&,const std::vector<ConnectivityCurve>&,double)const{}

// ── ConnectivityGenerator ─────────────────────────────────────
ConnectivityGenerator::ConnectivityGenerator(const LBFGSConfig&cfg):solver_(cfg){}

std::vector<SiblingCurve> ConnectivityGenerator::buildSiblings(
    const ConnId&id,const std::unordered_map<ConnId,BezierCurve>&done,const ClusterOrderSolver&cs)const
{
    std::vector<SiblingCurve>sibs;
    for(auto&[cid,curve]:done){if(cid==id)continue;
        SiblingCurve s; s.curve=curve;
        auto ex=cs.exemptionOf(id,cid);
        // If ClusterOrderSolver has no pair record (different entry lane),
        // it means these curves are from different road arms — structural
        // cross. Mark as exempt to avoid spurious cluster penalty.
        if (ex == CrossExemption::None) {
            // Check if there is actually a cluster pair between id and cid.
            // If not, this is a cross-traffic structural intersection → exempt.
            s.exempt_a1 = true;
        } else {
            s.exempt_a1=(ex==CrossExemption::StructuralCross);
        }
        s.exempt_a2_radius=(ex==CrossExemption::ObstacleCross)?1.5:0.0;
        sibs.push_back(std::move(s));}
    return sibs;
}

BezierCurve ConnectivityGenerator::postProcess(
    const BezierCurve&c,const SDFField&sdf,const Polygon2d&fence,double kmax,
    const Vec2d&t0_orig,const Vec2d&t1_orig,bool skip_elastic_band,
    const Vec2d* p0_exact,const Vec2d* p1_exact)
{
    // RC-3 FIX: after adaptiveRefine splits, warm-start re-optimise
    auto refined=adaptiveRefine(c,sdf,kmax);
    BezierCurve cur=refined.curve;
    if(refined.was_split){
        // Re-optimise the expanded curve with G1 tangent constraints
        PenaltyCost cost2;
        cost2.proto          = cur;
        cost2.sdf            = &sdf;
        cost2.fence          = fence;
        cost2.start_tan_dir  = t0_orig.norm()>1e-8 ? t0_orig.normalized() : c.startTan();
        cost2.end_tan_dir    = t1_orig.norm()>1e-8 ? t1_orig.normalized() : c.endTan();
        cost2.obstacle_clearance = 0.0;
        cost2.full_param_mode    = (cur.numSegments() > 1);
        cost2.buildCache();
        // Short warm-start: fewer outer iters, already near solution
        cur = optimiseCurve(cost2, solver_, cur, /*outer_iters=*/2);
    }

    // Derive exact endpoints and tangents from original lane geometry args
    Vec2d st = t0_orig.norm()>1e-8 ? t0_orig.normalized() : cur.startTan();
    Vec2d et = t1_orig.norm()>1e-8 ? t1_orig.normalized() : cur.endTan();

    // Use exact lane endpoint positions if provided; otherwise fall back to
    // the current curve endpoints.
    Vec2d ep0 = p0_exact ? *p0_exact : cur.startPt();
    Vec2d ep1 = p1_exact ? *p1_exact : cur.endPt();

    // For U-turns (skip_elastic_band=true) or any high-curvature curve where
    // elasticBandSmooth would produce oscillations: skip the band-smooth path
    // and just hard-pin endpoints + enforce G1 tangents analytically.
    if (skip_elastic_band) {
        if (!cur.segs.empty()) {
            // Hard-pin endpoints to exact lane geometry positions
            cur.segs.front().ctrl[0] = ep0;
            cur.segs.back() .ctrl[3] = ep1;
            // Re-enforce G1 tangent directions at endpoints.
            // Use a meaningful handle length based on the distance to the next
            // interior control point so the curve shape is preserved.
            {
                Vec2d& c1 = cur.segs.front().ctrl[1];
                Vec2d c2  = cur.segs.front().ctrl[2];
                // Compute handle length as projection of existing c1 onto st
                // but ensure it is at least seg_len * 0.1 for numerical stability
                double seg_len = (cur.segs.front().ctrl[3] - ep0).norm();
                double lam = (c1 - ep0).dot(st);
                lam = std::max(lam, seg_len * 0.1);
                c1 = ep0 + lam * st;
            }
            {
                Vec2d& c2 = cur.segs.back().ctrl[2];
                double seg_len = (ep1 - cur.segs.back().ctrl[0]).norm();
                double mu = (ep1 - c2).dot(et);
                mu = std::max(mu, seg_len * 0.1);
                c2 = ep1 - mu * et;
            }
        }
        return cur;
    }

    double arc = cur.arcLength();
    int n_samp  = std::max(40, (int)(arc / 0.15));
    auto pts = cur.sampleByArcLength(n_samp);
    if((int)pts.size()>=3){
        // BUG-1b FIX: hard-pin endpoints before smoothing
        pts.front() = ep0;
        pts.back()  = ep1;

        double max_k   = cur.maxCurvature(20);
        double move_st = std::min(0.15, std::max(0.02, max_k * 0.3));
        auto sm = elasticBandSmooth(pts,sdf,fence,kmax,move_st,80,0.1);

        // BUG-1c FIX: hard-pin endpoints after smoothing too
        // (elasticBand doesn't move i=0 and i=N-1, but float arithmetic
        //  may cause tiny drift; this guarantees exact match)
        sm.front() = ep0;
        sm.back()  = ep1;

        cur = rebuildFromSmoothedPts(sm, st, et);

        // BUG-1 FINAL: after rebuild, ensure ctrl[0] and ctrl[3] are exactly
        // at ep0 and ep1 (catmull-rom may drift by floating point)
        if(!cur.segs.empty()){
            cur.segs.front().ctrl[0] = ep0;
            cur.segs.back() .ctrl[3] = ep1;
        }
    }
    return cur;
}

void ConnectivityGenerator::validate(ConnectivityCurve&cc,const IntersectionInput&input,const SDFField&sdf)const{
    if(!cc.curve)return;
    auto&c=*cc.curve;
    double ms=minSDFAlongCurve(c,sdf,20);
    cc.violation.max_obstacle_penetration=std::max(0.0,-ms);
    if(!input.area.coarse_area.outer.empty()){
        double ov=0;
        for(auto&pt:c.sampleByArcLength(30))
            if(!polygonContains(input.area.coarse_area,pt))
                ov=std::max(ov,pointToPolygonDist(pt,input.area.coarse_area));
        cc.violation.max_fence_overflow=ov;
    }
    if(cc.violation.max_obstacle_penetration>0.05){
        cc.status=CurveStatus::Degraded;
        cc.violation.reason="Residual obstacle penetration after optimisation";
    } else if(!cc.violation.exempt_crosses.empty()){
        cc.status=CurveStatus::WarnA2;
    } else cc.status=CurveStatus::OK;
}

ConnectivityCurve ConnectivityGenerator::generateOne(
    const Connectivity&conn,const IntersectionInput&input,
    const SDFField&sdf,const SDFField&sdf_coarse,
    const std::vector<SiblingCurve>&siblings)
{
    ConnectivityCurve cc;
    cc.id=conn.id;cc.entry_lane_id=conn.entry_lane_id;
    cc.exit_lane_id=conn.exit_lane_id;cc.turn_type=conn.turn_type;

    auto[p0,t0]=input.entryPtDir(conn.entry_lane_id);
    auto[p1,t1]=input.exitPtDir(conn.exit_lane_id);

    double lw=3.5;
    if(auto*l=input.findLane(conn.entry_lane_id))lw=l->width;

    auto pre=preCheck(sdf_coarse,input.area.coarse_area,p0,p1,lw,input.boundaries);
    if(pre.type==ViolationInfo::InfeasibilityType::TopologicalBlock)
        return makeFallbackCurve(pre,conn,p0,p1);

    // Build sibling polylines from already-generated curves in this cluster.
    // Only include same-entry siblings (curves from the same entry lane) for the
    // initial curve shape selection.  Different-entry curves are cross-traffic
    // (structural intersection) and must not constrain initial curve shape.
    std::vector<std::vector<Vec2d>> sib_polys;
    for (auto& sib : siblings) {
        if (!sib.exempt_a1 && (sib.curve.startPt() - p0).norm() < 0.5) {
            sib_polys.push_back(sib.curve.sampleByArcLength(20));
        }
    }
    BezierCurve initial=buildInitialCurve(p0,t0,p1,t1,sdf,input.area.coarse_area,sib_polys);

    PenaltyCost cost;
    cost.proto=initial;cost.sdf=&sdf;
    cost.boundaries=input.boundaries;cost.fence=input.area.coarse_area;
    // Mark same-entry siblings as structurally exempt (they share an entry point
    // and their initial diverging overlap is geometrically expected / legal).
    // This prevents the cluster penalty from pushing turn curves to the outer
    // side just to avoid the overlap with the straight-through curve from the
    // same entry lane.
    auto modified_siblings = siblings;
    for (auto& sib : modified_siblings) {
        if (!sib.exempt_a1 && (sib.curve.startPt() - p0).norm() < 0.5) {
            sib.exempt_a1 = true;  // same-entry diverge: structural cross
        }
    }
    cost.siblings=modified_siblings;
    cost.crosswalks=input.crosswalks;
    // obstacle_clearance: match the penetration-only strategy used in
    // buildInitialCurve (INIT_CLEARANCE=0).  The SDF field has obstacle_buffer=0.4m
    // baked in, so SDF=0 means "touching the buffered obstacle surface" which
    // already represents being 0.4m from the raw obstacle edge.
    // Using clearance=0 means the optimizer only penalises actual penetration
    // (SDF<0), without trying to push curves away from obstacles that they
    // legitimately pass near (e.g. a straight lane 1.25m from an obstacle).
    // The 0.4m buffer already provides a meaningful safety margin.
    cost.obstacle_clearance=0.0;
    // G1 hard-constraint: pin endpoint tangent directions from Lane geometry
    cost.start_tan_dir=t0.norm()>1e-8?t0.normalized():Vec2d(1,0);
    cost.end_tan_dir  =t1.norm()>1e-8?t1.normalized():Vec2d(1,0);
    // Level-2: if initial curve has >1 segment (geometric bypass / multi-waypoint),
    // activate full-param mode so join points are also optimised.
    cost.full_param_mode = (initial.numSegments() > 1);

    BezierCurve opt=optimiseCurve(cost,solver_,initial,4);
    bool is_uturn = (conn.turn_type==TurnType::UTurnLeft||conn.turn_type==TurnType::UTurnRight);
    // Skip elasticBandSmooth for:
    // 1. U-turns (high curvature, elasticBand produces oscillations)
    // 2. Multi-segment curves (Level-1 arch, Level-2): elasticBand cannot preserve
    //    G1 at intermediate join points or at the endpoints.
    // 3. Single-segment curves where the initial Bezier arc already clears
    //    obstacles (SDF ≥ clearance): elasticBand would distort the G1.
    // In practice, since our optimizer enforces SDF clearance, elasticBand is
    // rarely needed and often harmful to G1 continuity.
    // We disable it globally and rely on the optimizer alone for smoothness.
    bool skip_band = true;  // always skip: preserve G1, rely on optimizer
    // Pass exact lane endpoint positions to guarantee G1 continuity regardless
    // of any optimiser drift during the multi-segment optimisation.
    BezierCurve final_c=postProcess(opt,sdf,input.area.coarse_area,0.25,t0,t1,skip_band,&p0,&p1);
    cc.curve=final_c;
    // Print final curve for QGIS verification
    std::cout << toWKT(genVectorline({p0[0],p0[1],0}, {t0[0],t0[1],0}, 1.0), conn.entry_lane_id) << std::endl;
    std::cout << toWKT(genVectorline({p1[0],p1[1],0}, {t1[0],t1[1],0}, 1.0), conn.exit_lane_id) << std::endl;
    std::cout << toWKT(toArray(final_c.sample(50)), cc.id) << std::endl;
    validate(cc,input,sdf);
    if(pre.narrow_passage&&cc.status==CurveStatus::OK)cc.status=CurveStatus::WarnA2;
    return cc;
}

std::vector<ConnectivityCurve> ConnectivityGenerator::generate(
    const IntersectionInput&input,SDFField&sdf,double*out_ms)
{
    auto t0=std::chrono::steady_clock::now();

    SDFField sdf_coarse;
    auto roi=input.area.coarse_area.empty()?BoundingBox2d{}:input.area.coarse_area.bbox();
    if(roi.width()<1){
        for(auto&l:input.lanes)
            for(auto&p:l.geometry.points)roi.expand(p);
        //if(roi.empty())for(auto&g:input.lane_groups)roi.expand(g.ref_point);
        roi.min_pt-=Vec2d(20,20);roi.max_pt+=Vec2d(20,20);}
    sdf_coarse.build(roi,input.obstacles,0.5,0.4);

    cluster_solver_.build(input.connectivities,input.lanes,input.lane_groups);

    GlobalCoordinator coord;
    coord.build(input.connectivities,input);

    std::vector<ConnectivityCurve>results;
    results.reserve(input.connectivities.size());

    std::unordered_map<ConnId,const Connectivity*>cmap;
    for(auto&c:input.connectivities)cmap[c.id]=&c;

    std::unordered_map<ConnId,BezierCurve>done;

    for(auto&group:coord.groups()){
        for(auto&cid:group.conn_ids){
            auto*conn=cmap[cid]; if(!conn)continue;
            auto sibs=buildSiblings(cid,done,cluster_solver_);
            auto cc=generateOne(*conn,input,sdf,sdf_coarse,sibs);
            if(cc.curve)done[cid]=*cc.curve;
            results.push_back(std::move(cc));
        }
        cluster_solver_.checkAndMarkA2(done,sdf,1.5);
    }

    auto t1=std::chrono::steady_clock::now();
    if(out_ms)*out_ms=std::chrono::duration<double,std::milli>(t1-t0).count();
    return results;
}
