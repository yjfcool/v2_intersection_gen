#pragma once
#include "types.h"
class FineAreaGenerator{
public:
    Polygon2d generate(const IntersectionInput&,const std::vector<ConnectivityCurve>&);
private:
    std::vector<LineString2d> sortEdges(const std::vector<LineString2d>&,double snap=0.1)const;
    std::vector<Vec2d> collectBoundaryPts(const std::vector<LineString2d>&)const;
    static Polygon2d convexHull(const std::vector<Vec2d>&);
    static bool isSimple(const Polygon2d&);
};
