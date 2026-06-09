#include "connectivity_generator.h"
#include "curve/hermite_init.h"
#include "curve/curve_utils.h"
#include "constraints/fence_check.h"
#include "utils.h"
#include "optimizer/sdf_field.h"
#include "constraints/infeasibility_detector.h"
#include "utils/quadtree.h"
#include <chrono>
#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <memory>

// ── IntersectionInput helpers ─────────────────────────────────────────────────
const bool IntersectionInput::IsEntryLane(const LaneId& id) const {
    for (auto& lg : lane_groups)
        if (std::find(lg.lanes.begin(), lg.lanes.end(), id) != lg.lanes.end())
            return lg.role == GroupRole::Entry;
    return false;
}

const bool IntersectionInput::IsEntryLaneEdge(const LaneEdgeId& id) const {
    for (auto& lg : lane_groups)
        if (std::find(lg.boundaries.begin(), lg.boundaries.end(), id) != lg.boundaries.end())
            return lg.role == GroupRole::Entry;
    for (auto& lg : lane_groups)
        for (auto& lid : lg.lanes) {
            auto l = findLane(lid);
            if (l && (l->left_edge_id == id || l->right_edge_id == id))
                return lg.role == GroupRole::Entry;
        }
    return false;
}

const Lane* IntersectionInput::findLane(const LaneId& id) const {
    for (auto& l : lanes) if (l.id == id)return &l;
    return nullptr;
}

const LaneGroup* IntersectionInput::findGroup(const LaneGroupId& id) const {
    for (auto& g : lane_groups) if (g.id == id)return &g;
    return nullptr;
}

const LaneEdge* IntersectionInput::findEdge(const LaneEdgeId& id) const {
    for (auto& e : lane_edges) if (e.id == id)return &e;
    return nullptr;
}

bool IntersectionInput::laneGroupExists(const LaneGroupId& id) const {
    return findGroup(id) != nullptr;
}

std::pair<Vec2d, Vec2d> IntersectionInput::entryPtDir(const LaneId& lid) const {
    auto* l = findLane(lid);
    if (!l || l->geometry.points.empty()) {
        std::cout << "[WARN] entrylane:" << lid << " no geometry!\n";
        return {Vec2d(0, 0), Vec2d(1, 0)};
    }
    return {entryLinePoint(l->geometry.points), entryLineTangent(l->geometry.points)};
}

std::pair<Vec2d, Vec2d> IntersectionInput::exitPtDir(const LaneId& lid) const {
    auto* l = findLane(lid);
    if (!l || l->geometry.points.empty()) {
        std::cout << "[WARN] exitlane:" << lid << " no geometry!\n";
        return {Vec2d(10, 0), Vec2d(1, 0)};
    }
    return {exitLinePoint(l->geometry.points), exitLineTangent(l->geometry.points)};
}

// ── GlobalCoordinator ─────────────────────────────────────────────────────────
// Generation order: straight first (anchor), then left/right turns,
// then U-turns last (they are most constrained and need all others as siblings).
static int globalTurnPriority(ConnTurnType t) {
    switch (t) {
    case ConnTurnType::Straight: return 0;
    case ConnTurnType::TurnLeft: return 1;
    case ConnTurnType::TurnRight: return 1;
    case ConnTurnType::UTurnLeft: return 2;
    case ConnTurnType::UTurnRight: return 2;
    default: return 0;
    }
}

void GlobalCoordinator::build(
    const std::vector<Connectivity>& conns, const IntersectionInput& inp) {
    std::map<int, OptGroup> pm;
    for (auto& c : conns) {
        int p = globalTurnPriority(c.turn_type);
        pm[p].conn_ids.push_back(c.id);
        pm[p].priority = p;
    }
    groups_.clear();
    for (auto& kv : pm) {
        groups_.push_back(kv.second);
    }
}

void GlobalCoordinator::addSoftObstacles(SDFField&,
    const std::vector<ConnectivityCurve>&, double) const {}

// ── ConnectivityGenerator ─────────────────────────────────────────────────────

void ConnectivityGenerator::initializePersistentQuadTree(const IntersectionInput& input) {
    if (!quad_tree_) {
        // Determine overall bounds from the intersection area
        BoundingBox2d bounds;
        if (!input.area.geometry.outer.empty()) {
            bounds = input.area.geometry.bbox();
        } else {
            // If no area defined, calculate from lane geometries
            for (const auto& lane : input.lanes) {
                for (const auto& pt : lane.geometry.points) {
                    bounds.expand(pt);
                }
            }
            // Add buffer
            bounds.min_pt -= Vec2d(50.0, 50.0);
            bounds.max_pt += Vec2d(50.0, 50.0);
        }

        // Initialize the persistent quadtree with larger capacity
        quad_tree_ = make_unique_cpp11<QuadTree>(bounds, 16, 8); // increased capacity and depth for better performance
    }
}

void ConnectivityGenerator::addCurveToPersistentQuadTree(const ConnId& id, const BezierCurve& curve) {
    if (quad_tree_) {
        quad_tree_->insert(curve, id);
    }
}

void ConnectivityGenerator::clearPersistentSpatialIndex() {
    quad_tree_.reset();
}

ConnectivityGenerator::ConnectivityGenerator(const LBFGSConfig& cfg) : solver_(cfg), quad_tree_(nullptr) {}

std::vector<SiblingCurve> ConnectivityGenerator::buildSiblingsOptimized(
    const ConnId& id, const std::unordered_map<ConnId, BezierCurve>& done,
    const ClusterOrderSolver& cs, const IntersectionInput& input) const {

    std::vector<SiblingCurve> sibs;

    // Find the current connection to get its entry lane for spatial query
    const Connectivity* current_conn = nullptr;
    for (const auto& conn : input.connectivities) {
        if (conn.id == id) {
            current_conn = &conn;
            break;
        }
    }

    if (!current_conn || done.empty()) {
        return buildSiblings(id, done, cs); // fallback to original if needed
    }

    // Determine query bounds based on current connection
    BoundingBox2d query_bounds = getCurveBoundsForQuery(*current_conn, input);

    // Use the persistent quadtree instance from the class if available
    if (this->quad_tree_) {
        auto nearby_curves = this->quad_tree_->queryRange(query_bounds);

        for (const auto& kv : nearby_curves) {
            auto& curve = kv.first;
            auto& conn_id = kv.second;

            // Skip the curve with the same ID
            if (conn_id == id) continue;

            SiblingCurve s;
            s.curve = curve;

            auto ex = cs.exemptionOf(id, conn_id);
            if (ex == CrossExemption::None) {
                s.exempt_a1 = true; // cross-traffic, different entry & exit → structural cross
            } else {
                s.exempt_a1 = (ex == CrossExemption::StructuralCross);
            }

            s.exempt_a2_radius = (ex == CrossExemption::ObstacleCross) ? 1.5 : 0.0;
            sibs.push_back(std::move(s));
        }
    } else {
        // Fallback to temporary quadtree construction if persistent one is not available
        QuadTree temp_quadtree(query_bounds, 8, 6); // capacity=8, max depth=6

        for (const auto& kv : done) {
            auto& conn_id = kv.first;
            auto& curve = kv.second;
            if (conn_id == id) continue;
            temp_quadtree.insert(curve, conn_id);
        }

        // Query for potentially overlapping curves
        auto nearby_curves = temp_quadtree.queryRange(query_bounds);

        for (const auto& kv : nearby_curves) {
            auto& curve = kv.first;
            auto& conn_id = kv.second;
            SiblingCurve s;
            s.curve = curve;

            auto ex = cs.exemptionOf(id, conn_id);
            if (ex == CrossExemption::None) {
                s.exempt_a1 = true; // cross-traffic, different entry & exit → structural cross
            } else {
                s.exempt_a1 = (ex == CrossExemption::StructuralCross);
            }

            s.exempt_a2_radius = (ex == CrossExemption::ObstacleCross) ? 1.5 : 0.0;
            sibs.push_back(std::move(s));
        }
    }

    return sibs;
}

std::vector<SiblingCurve> ConnectivityGenerator::buildSiblings(
    const ConnId& id, const std::unordered_map<ConnId, BezierCurve>& done, const ClusterOrderSolver& cs) const {
    std::vector<SiblingCurve> sibs;
    for (auto& kv : done) {
        auto& cid = kv.first;
        auto& curve = kv.second;
        if (cid == id)
            continue;

        SiblingCurve s;
        s.curve = curve;

        auto ex = cs.exemptionOf(id, cid);
        // ClusterOrderSolver now registers pairs for BOTH same-entry AND same-exit
        // clusters. If exemptionOf returns None, there is genuinely no cluster
        // relationship (different entry AND different exit lanes) → structural
        // cross-traffic from different road arms → exempt from cluster penalty.
        if (ex == CrossExemption::None) {
            s.exempt_a1 = true; // cross-traffic, different entry & exit → structural cross
        } else {
            s.exempt_a1 = (ex == CrossExemption::StructuralCross);
        }

        s.exempt_a2_radius = (ex == CrossExemption::ObstacleCross) ? 1.5 : 0.0;
        sibs.push_back(std::move(s));
    }
    return sibs;
}

BoundingBox2d ConnectivityGenerator::getCurveBoundsForQuery(const Connectivity& conn, const IntersectionInput& input) const {
    // Get entry and exit points for the connection
    auto _entry = input.entryPtDir(conn.entry_lane_id);
    Vec2d p0 = _entry.first;
    Vec2d t0 = _entry.second;
    auto _exit = input.exitPtDir(conn.exit_lane_id);
    Vec2d p1 = _exit.first;

    // Create a bounding box centered around the connection with a reasonable buffer
    BoundingBox2d bounds;
    bounds.min_pt = Vec2d(std::min(p0.x(), p1.x()), std::min(p0.y(), p1.y()));
    bounds.max_pt = Vec2d(std::max(p0.x(), p1.x()), std::max(p0.y(), p1.y()));

    // Add buffer to account for curve expansion
    bounds.min_pt -= Vec2d(10.0, 10.0);
    bounds.max_pt += Vec2d(10.0, 10.0);

    return bounds;
}


BezierCurve ConnectivityGenerator::postProcess(
    const BezierCurve& c, const SDFField& sdf, const Polygon2d& fence, double kmax,
    const Vec2d& t0_orig, const Vec2d& t1_orig, bool skip_elastic_band,
    const Vec2d* p0_exact, const Vec2d* p1_exact) {
    auto refined = adaptiveRefine(c, sdf, kmax);
    BezierCurve cur = refined.curve;
    if (refined.was_split) {
        PenaltyCost cost2;
        cost2.proto = cur;
        cost2.sdf = &sdf;
        cost2.fence = fence;
        cost2.start_tan_dir = t0_orig.norm() > 1e-8 ? t0_orig.normalized() : c.startTan();
        cost2.end_tan_dir = t1_orig.norm() > 1e-8 ? t1_orig.normalized() : c.endTan();
        cost2.obstacle_clearance = 0.0;
        cost2.full_param_mode = (cur.numSegments() > 1);
        cost2.buildCache();
        cur = optimiseCurve(cost2, solver_, cur, 2);
    }

    Vec2d st = t0_orig.norm() > 1e-8 ? t0_orig.normalized() : cur.startTan();
    Vec2d et = t1_orig.norm() > 1e-8 ? t1_orig.normalized() : cur.endTan();
    Vec2d ep0 = p0_exact ? *p0_exact : cur.startPt();
    Vec2d ep1 = p1_exact ? *p1_exact : cur.endPt();

    if (skip_elastic_band) {
        if (!cur.segs.empty()) {
            cur.segs.front().ctrl[0] = ep0;
            cur.segs.back().ctrl[3] = ep1;
            Vec2d& c1 = cur.segs.front().ctrl[1];
            double sl = (cur.segs.front().ctrl[3] - ep0).norm();
            double lam = std::max((c1 - ep0).dot(st), sl * 0.1);
            c1 = ep0 + lam * st;
            Vec2d& c2 = cur.segs.back().ctrl[2];
            double el = (ep1 - cur.segs.back().ctrl[0]).norm();
            double mu = std::max((ep1 - c2).dot(et), el * 0.1);
            c2 = ep1 - mu * et;
        }
        return cur;
    }

    double arc = cur.arcLength();
    int n_samp = std::max(40, (int)(arc / 0.15));
    auto pts = cur.sampleByArcLength(n_samp);
    if ((int)pts.size() >= 3) {
        pts.front() = ep0;
        pts.back() = ep1;

        double max_k = cur.maxCurvature(20);
        double move_st = std::min(0.15, std::max(0.02, max_k * 0.3));
        // Using adaptive version of elastic band smoothing
        auto sm = elasticBandSmoothAdaptive(pts, sdf, fence, kmax, move_st, 80, 0.1);
        sm.front() = ep0;
        sm.back() = ep1;

        cur = rebuildFromSmoothedPts(sm, st, et);
        if (!cur.segs.empty()) {
            cur.segs.front().ctrl[0] = ep0;
            cur.segs.back().ctrl[3] = ep1;
        }
    }
    return cur;
}

void ConnectivityGenerator::validate(
    ConnectivityCurve& cc, const IntersectionInput& input, const SDFField& sdf) const {
    if (!cc.curve) return;
    auto& c = *cc.curve;
    double ms = minSDFAlongCurve(c, sdf, 20);
    cc.violation.max_obstacle_penetration = std::max(0.0, -ms);
    if (!input.area.geometry.outer.empty()) {
        double ov = 0;
        for (auto& pt : c.sampleByArcLength(30))
            if (!polygonContains(input.area.geometry, pt))
                ov = std::max(ov, pointToPolygonDist(pt, input.area.geometry));
        cc.violation.max_fence_overflow = ov;
    }
    if (cc.violation.max_obstacle_penetration > 0.05)
        cc.status = CurveStatus::Degraded;
    else if (!cc.violation.exempt_crosses.empty())
        cc.status = CurveStatus::WarnA2;
    else
        cc.status = CurveStatus::OK;
}

ConnectivityCurve ConnectivityGenerator::generateOne(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf,
    const SDFField& sdf_coarse, const std::vector<SiblingCurve>& siblings) {
    ConnectivityCurve cc;
    cc.id = conn.id;
    cc.entry_lane_id = conn.entry_lane_id;
    cc.exit_lane_id = conn.exit_lane_id;
    cc.turn_type = conn.turn_type;

    auto _entry = input.entryPtDir(conn.entry_lane_id);
    Vec2d p0 = _entry.first;
    Vec2d t0 = _entry.second;
    auto _exit = input.exitPtDir(conn.exit_lane_id);
    Vec2d p1 = _exit.first;
    Vec2d t1 = _exit.second;

    double lw = 3.5;
    if (auto* l = input.findLane(conn.entry_lane_id))
        lw = l->width;

    auto pre = preCheck(sdf_coarse, input.area.geometry, p0, p1, lw, input.boundaries);
    if (pre.type == ViolationInfo::InfeasibilityType::TopologicalBlock)
        return makeFallbackCurve(pre, conn, p0, p1);

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
    BezierCurve initial = buildInitialCurve(p0, t0, p1, t1, sdf, input.area.geometry, sib_polys);

    // ── Optimise ──────────────────────────────────────────────────────────────
    PenaltyCost cost;
    cost.proto = initial;
    cost.sdf = &sdf;
    cost.boundaries = input.boundaries;
    cost.fence = input.area.geometry;
    cost.siblings = siblings;
    cost.crosswalks = input.crosswalks;
    cost.obstacle_clearance = 0.0;
    cost.start_tan_dir = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    cost.end_tan_dir = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    cost.full_param_mode = (initial.numSegments() > 1);

    BezierCurve opt = optimiseCurveWithEarlyStopping(cost, solver_, initial, /*outer_iters=*/4);

    bool skip_band = true; // always skip: preserve G1, rely on optimizer
    BezierCurve final_c = postProcess(opt, sdf, input.area.geometry, 0.25, t0, t1, skip_band, &p0, &p1);
    cc.curve = std::make_shared<BezierCurve>(final_c);
    validate(cc, input, sdf);
    if (pre.narrow_passage && cc.status == CurveStatus::OK)
        cc.status = CurveStatus::WarnA2;
    return cc;
}

std::vector<ConnectivityCurve> ConnectivityGenerator::generate(
    const IntersectionInput& input, SDFField& sdf, double* out_ms) {
    auto t0 = std::chrono::steady_clock::now();

    SDFField sdf_coarse;
    auto roi = input.area.geometry.empty() ? BoundingBox2d{} : input.area.geometry.bbox();
    if (roi.width() < 1) {
        for (auto& l : input.lanes) {
            for (auto& p : l.geometry.points)
                roi.expand(p);
        }
        roi.min_pt -= Vec2d(20, 20);
        roi.max_pt += Vec2d(20, 20);
    }
    sdf_coarse.build(roi, input.obstacles, 0.5, 0.4);

    // Build group-based cluster solver
    cluster_solver_.build(input.connectivities, input.lanes, input.lane_groups);

    // Build generation order
    GlobalCoordinator coord;
    coord.build(input.connectivities, input);

    // Initialize persistent spatial index for efficient sibling queries
    initializePersistentQuadTree(input);

    std::vector<ConnectivityCurve> results;
    results.reserve(input.connectivities.size());
    std::unordered_map<ConnId, const Connectivity*> cmap;
    for (auto& c : input.connectivities) {
        cmap[c.id] = &c;
    }
    std::unordered_map<ConnId, BezierCurve> done;
    for (auto& group : coord.groups()) {
        for (auto& cid : group.conn_ids) {
            auto* conn = cmap[cid];
            if (!conn) continue;

            // Use optimized version with spatial indexing
            auto sibs = buildSiblingsOptimized(cid, done, cluster_solver_, input);

            auto cc = generateOne(*conn, input, sdf, sdf_coarse, sibs);
            if (cc.curve) {
                done[cid] = *cc.curve;

                // Add the newly generated curve to the persistent spatial index
                addCurveToPersistentQuadTree(cid, *cc.curve);
            }

            results.push_back(std::move(cc));
        }
        // After each priority group: mark obstacle-adjacent crossings as soft
        cluster_solver_.checkAndMarkA2(done, sdf, 1.5);
    }

    // Clean up the persistent spatial index after processing
    clearPersistentSpatialIndex();

    auto t1 = std::chrono::steady_clock::now();
    if (out_ms) *out_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return results;
}
