#pragma once
#include "types.h"
#include "curve/bezier.h"
#include "optimizer/sdf_field.h"
class LaneEdgeGenerator{
public:
    std::vector<ConnectivityLaneEdge> generate(const IntersectionInput&,const std::vector<ConnectivityCurve>&,const SDFField&);
private:
    const BezierCurve* curveForLane(const LaneId&,const std::vector<ConnectivityCurve>&)const;
    ConnectivityLaneEdge makeSharedEdge(const LaneEdgeId&,const BezierCurve&,const BezierCurve&,const Vec2d&,const Vec2d&,bool);
    ConnectivityLaneEdge makeOffsetEdge(const LaneEdgeId&,const BezierCurve&,double,const Vec2d&,const Vec2d&);
};
