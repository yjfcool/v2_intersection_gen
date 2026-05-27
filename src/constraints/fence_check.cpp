#include "fence_check.h"
#include "curve/bezier.h"
#include "curve/curve_utils.h"
#include <cmath>
#include <algorithm>

// Polygon2d::contains - definition for the method declared in types.h
bool Polygon2d::contains(const Vec2d&p)const{return polygonContains(*this,p);}

bool polygonContains(const Polygon2d&poly,const Vec2d&pt){
    auto&ring=poly.outer;int n=(int)ring.size();if(n<3)return false;
    bool inside=false;
    for(int i=0,j=n-1;i<n;j=i++){
        double xi=ring[i].x(),yi=ring[i].y(),xj=ring[j].x(),yj=ring[j].y();
        if(((yi>pt.y())!=(yj>pt.y()))&&(pt.x()<(xj-xi)*(pt.y()-yi)/(yj-yi+1e-20)+xi))
            inside=!inside;}
    if(inside)for(auto&h:poly.holes){Polygon2d hp;hp.outer=h;if(polygonContains(hp,pt))return false;}
    return inside;}
double pointToPolygonDist(const Vec2d&pt,const Polygon2d&poly){
    auto&ring=poly.outer;int n=(int)ring.size();if(n<2)return 1e18;
    double m=1e18;
    for(int i=0;i<n;++i){
        Vec2d a=ring[i],b=ring[(i+1)%n],ab=b-a,ap=pt-a;
        double t=std::max(0.0,std::min(1.0,ab.dot(ap)/std::max(1e-20,ab.squaredNorm())));
        m=std::min(m,(pt-(a+t*ab)).norm());}
    return m;}
bool curveIntersectsBoundary(const BezierCurve&c,const std::vector<Boundary>&bds,int sps){
    std::vector<Vec2d>cp;
    for(auto&seg:c.segs)for(int i=0;i<=sps;++i)cp.push_back(seg.evaluate((double)i/sps));
    for(auto&bnd:bds)for(int bi=0;bi+1<(int)bnd.geometry.points.size();++bi)
        for(int ci=0;ci+1<(int)cp.size();++ci)
            if(segmentsIntersect(cp[ci],cp[ci+1],bnd.geometry.points[bi],bnd.geometry.points[bi+1]))
                return true;
    return false;}
bool curveInsideFence(const BezierCurve&c,const Polygon2d&fence,int sps){
    if(fence.outer.empty())return true;
    for(auto&seg:c.segs)for(int i=0;i<=sps;++i)
        if(!polygonContains(fence,seg.evaluate((double)i/sps)))return false;
    return true;}
