#include "lane_edge_generator.h"
#include "curve/curve_utils.h"
#include "utils.h"

static std::vector<Vec2d> offsetPL(const std::vector<Vec2d>&pts,double off){
    std::vector<Vec2d>out; out.reserve(pts.size()); int n=(int)pts.size();
    for(int i=0;i<n;++i){
        Vec2d tan;
        if(i==0)tan=(pts[1]-pts[0]).normalized();
        else if(i==n-1)tan=(pts[n-1]-pts[n-2]).normalized();
        else tan=(pts[i+1]-pts[i-1]).normalized();
        Vec2d nm{-tan[1],tan[0]};
        out.push_back(pts[i]+off*nm);
    }
    return out;
}

const BezierCurve* LaneEdgeGenerator::curveForLane(const LaneId&lid,const std::vector<ConnectivityCurve>&ccs)const{
    for(auto&cc:ccs)if(cc.entry_lane_id==lid&&cc.curve)return&(*cc.curve);
    return nullptr;
}

ConnectivityLaneEdge LaneEdgeGenerator::makeSharedEdge(const LaneEdgeId&eid,const BezierCurve&a,const BezierCurve&b,
    const Vec2d&et,const Vec2d&xt,bool shared)
{
    auto mid=midlineSampleByArcLength(a,b,20);
    auto ec=fitBezierWithEndTangents(mid,et,xt);
    ConnectivityLaneEdge e; e.id=eid; e.is_shared=shared;
    e.geometry.points=ec.sampleByArcLength(30);
    return e;
}

ConnectivityLaneEdge LaneEdgeGenerator::makeOffsetEdge(const LaneEdgeId&eid,const BezierCurve&c,double off,const Vec2d&,const Vec2d&){
    auto pts=c.sampleByArcLength(30);
    ConnectivityLaneEdge e; e.id=eid; e.is_shared=false;
    e.geometry.points=offsetPL(pts,off);
    return e;
}

std::vector<ConnectivityLaneEdge> LaneEdgeGenerator::generate(
    const IntersectionInput&input,const std::vector<ConnectivityCurve>&ccs,const SDFField&)
{
    std::vector<ConnectivityLaneEdge>result;
    for(auto&group:input.lane_groups){
        if(group.lanes.empty())continue;
        int N=(int)group.lanes.size();
        std::vector<const BezierCurve*>lc(N,nullptr);
        for(int i=0;i<N;++i)lc[i]=curveForLane(group.lanes[i],ccs);
        for(int b=0;b<=N;++b){
            LaneEdgeId eid=group.id+"_ie_"+std::to_string(b);
            const BezierCurve*lft=(b>0)?lc[b-1]:nullptr;
            const BezierCurve*rgt=(b<N)?lc[b]:nullptr;

            // Derive entry/exit tangents from Lane::geometry (not LaneGroup::direction)
            auto getLaneTangents=[&](int li) -> std::pair<Vec2d,Vec2d> {
                if(li<0||li>=N) return{{1,0},{1,0}};//{group.direction,group.direction};
                auto*ln=input.findLane(group.lanes[li]);
                if(!ln||ln->geometry.points.empty())
                    return{{1,0},{1,0}};//{group.direction,group.direction};
                if(group.role==GroupRole::Entry)
                    return{entryLineTangent(ln->geometry.points),entryLineTangent(ln->geometry.points)};
                else
                    return{exitLineTangent(ln->geometry.points),exitLineTangent(ln->geometry.points)};
            };
            // Use average of adjacent lane tangents for shared boundaries
            std::pair<Vec2d,Vec2d> _tl=getLaneTangents(b>0?b-1:b); Vec2d et_l=_tl.first;
            std::pair<Vec2d,Vec2d> _tr=getLaneTangents(b<N?b:b-1); Vec2d et_r=_tr.first;
            Vec2d et=(et_l+et_r);
            if(et.norm()<1e-8)et={1,0};//group.direction;
            else et.normalize();
            Vec2d xt=et;

            if(b>0&&b<N&&lft&&rgt){
                auto e=makeSharedEdge(eid,*lft,*rgt,et,xt,true);
                e.shared_by=std::make_shared<std::pair<LaneId, LaneId>>(group.lanes[b-1],group.lanes[b]);
                result.push_back(std::move(e));
            } else if(b==0&&rgt){
                double hw=1.75;if(auto*l=input.findLane(group.lanes[0]))hw=l->width*0.5;
                result.push_back(makeOffsetEdge(eid,*rgt,hw,et,xt));
            } else if(b==N&&lft){
                double hw=1.75;if(auto*l=input.findLane(group.lanes[N-1]))hw=l->width*0.5;
                result.push_back(makeOffsetEdge(eid,*lft,-hw,et,xt));
            }
        }
    }
    return result;
}
