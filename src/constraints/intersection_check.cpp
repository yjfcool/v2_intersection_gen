#include "intersection_check.h"
#include "curve/curve_utils.h"
bool bboxOverlap(const BezierCurve&a,const BezierCurve&b){return a.bbox().intersects(b.bbox());}
static void segPairCross(const BezierSegment&a,const BezierSegment&b,
    std::vector<Vec2d>&out,double tol,int depth=0)
{
    if(!a.bbox().intersects(b.bbox()))return;
    if(depth>10){out.push_back(a.evaluate(0.5));return;}
    if(a.arcLength(4)<tol&&b.arcLength(4)<tol){out.push_back(a.evaluate(0.5));return;}
    auto[al,ar]=a.splitAt(0.5);auto[bl,br]=b.splitAt(0.5);
    segPairCross(al,bl,out,tol,depth+1);segPairCross(al,br,out,tol,depth+1);
    segPairCross(ar,bl,out,tol,depth+1);segPairCross(ar,br,out,tol,depth+1);}
std::vector<Vec2d> curveCrossings(const BezierCurve&a,const BezierCurve&b,double tol){
    std::vector<Vec2d>raw;if(!bboxOverlap(a,b))return raw;
    for(auto&sa:a.segs)for(auto&sb:b.segs)segPairCross(sa,sb,raw,tol);
    std::vector<Vec2d>out;
    for(auto&p:raw){bool dup=false;
        for(auto&q:out)if((p-q).norm()<tol){dup=true;break;}
        if(!dup)out.push_back(p);}
    return out;}
bool curvesIntersectBusiness(const BezierCurve&a,const BezierCurve&b,double ep){
    auto cr=curveCrossings(a,b);
    for(auto&pt:cr)if(distToAllEndpoints(pt,a,b)>ep)return true;
    return false;}
