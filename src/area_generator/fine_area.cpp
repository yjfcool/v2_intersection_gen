#include "fine_area.h"
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include <algorithm>
#include <cmath>

std::vector<LineString2d> FineAreaGenerator::sortEdges(const std::vector<LineString2d>&edges,double snap)const{
    if(edges.empty())return{};
    std::vector<bool>used(edges.size(),false);
    std::vector<LineString2d>ordered;
    ordered.push_back(edges[0]);used[0]=true;
    for(int iter=0;iter<(int)edges.size()-1;++iter){
        Vec2d tail=ordered.back().points.back();
        double bd=1e18;int bi=-1;bool bf=false;
        for(int i=0;i<(int)edges.size();++i){
            if(used[i])continue;auto&pts=edges[i].points;if(pts.empty())continue;
            double dh=(pts.front()-tail).norm(),dt=(pts.back()-tail).norm();
            if(dh<bd){bd=dh;bi=i;bf=false;}
            if(dt<bd){bd=dt;bi=i;bf=true;}
        }
        if(bi<0||bd>snap*20)break;
        used[bi]=true;LineString2d ls=edges[bi];
        if(bf)std::reverse(ls.points.begin(),ls.points.end());
        ordered.push_back(ls);
    }
    for(int i=0;i<(int)edges.size();++i)if(!used[i])ordered.push_back(edges[i]);
    return ordered;
}

std::vector<Vec2d> FineAreaGenerator::collectBoundaryPts(const std::vector<LineString2d>&sorted)const{
    std::vector<Vec2d>pts;
    for(auto&ls:sorted)for(auto&p:ls.points){
        if(!pts.empty()&&(p-pts.back()).norm()<0.05)continue;
        pts.push_back(p);}
    if(!pts.empty()&&(pts.front()-pts.back()).norm()>0.1)pts.push_back(pts.front());
    return pts;
}

Polygon2d FineAreaGenerator::convexHull(const std::vector<Vec2d>&pts){
    if(pts.size()<3){Polygon2d p;p.outer=pts;return p;}
    std::vector<Vec2d>s=pts;
    auto pivot=*std::min_element(s.begin(),s.end(),[](const Vec2d&a,const Vec2d&b){
        return a[1]<b[1]||(a[1]==b[1]&&a[0]<b.x());});
    std::sort(s.begin(),s.end(),[&](const Vec2d&a,const Vec2d&b){
        double ax=a[0]-pivot[0],ay=a[1]-pivot[1],bx=b[0]-pivot[0],by=b[1]-pivot[1];
        double cr=ax*by-ay*bx;
        if(std::abs(cr)>1e-10)return cr>0;
        return ax*ax+ay*ay<bx*bx+by*by;});
    std::vector<Vec2d>hull;
    for(auto&p:s){
        while(hull.size()>=2){
            Vec2d a=hull[hull.size()-2],b=hull.back();
            if(cross2d(b-a,p-a)<=0)hull.pop_back();else break;}
        hull.push_back(p);}
    Polygon2d out;out.outer=hull;return out;
}

bool FineAreaGenerator::isSimple(const Polygon2d&poly){
    auto&ring=poly.outer;int n=(int)ring.size();if(n<4)return true;
    for(int i=0;i<n;++i)for(int j=i+2;j<n;++j){
        if(i==0&&j==n-1)continue;
        if(segmentsIntersect(ring[i],ring[(i+1)%n],ring[j],ring[(j+1)%n]))return false;}
    return true;
}

Polygon2d FineAreaGenerator::generate(const IntersectionInput&input,const std::vector<ConnectivityCurve>&ccs){
    std::vector<LineString2d>road_edges;
    for(auto&bnd:input.boundaries)if(bnd.type==BoundaryType::RoadEdge)road_edges.push_back(bnd.geometry);
    std::vector<Vec2d>endpoints;
    for(auto&cc:ccs){if(!cc.curve)continue;endpoints.push_back(cc.curve->startPt());endpoints.push_back(cc.curve->endPt());}
    if(road_edges.empty()){
        if(endpoints.empty())return input.area.coarse_area;
        auto pts=input.area.coarse_area.outer;pts.insert(pts.end(),endpoints.begin(),endpoints.end());
        return convexHull(pts);}
    auto sorted=sortEdges(road_edges,0.1);
    auto bpts=collectBoundaryPts(sorted);
    for(auto&ep:endpoints){
        double bd=1e18;int bi=0;
        for(int i=0;i+1<(int)bpts.size();++i){
            Vec2d ab=bpts[i+1]-bpts[i],ap=ep-bpts[i];
            double t=std::max(0.0,std::min(1.0,ab.dot(ap)/std::max(1e-12,ab.squaredNorm())));
            double d=(ep-(bpts[i]+t*ab)).norm();if(d<bd){bd=d;bi=i;}}
        if(bd<2.0)bpts.insert(bpts.begin()+bi+1,ep);
    }
    Polygon2d result;result.outer=bpts;
    if(result.outer.size()>1&&(result.outer.front()-result.outer.back()).norm()<0.05)result.outer.pop_back();
    if(!isSimple(result)){
        auto ap=bpts;ap.insert(ap.end(),endpoints.begin(),endpoints.end());
        result=convexHull(ap);}
    return result;
}
