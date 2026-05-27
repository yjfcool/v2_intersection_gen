#pragma once
#include "types.h"
bool polygonContains(const Polygon2d&,const Vec2d&);
double pointToPolygonDist(const Vec2d&,const Polygon2d&);
bool curveIntersectsBoundary(const BezierCurve&,const std::vector<Boundary>&,int sps=25);
bool curveInsideFence(const BezierCurve&,const Polygon2d&,int sps=25);
