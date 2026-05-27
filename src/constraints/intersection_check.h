#pragma once
#include "types.h"
std::vector<Vec2d> curveCrossings(const BezierCurve&,const BezierCurve&,double tol=0.01);
bool bboxOverlap(const BezierCurve&,const BezierCurve&);
bool curvesIntersectBusiness(const BezierCurve&,const BezierCurve&,double ep=0.01);
