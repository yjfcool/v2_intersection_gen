#include "sdf_field.h"
#include "constraints/fence_check.h"
#include <clipper2/clipper.h>
#include <cmath>
#include <algorithm>
#include <limits>

static Clipper2Lib::PathD toCP(const std::vector<Vec2d>&pts){
    Clipper2Lib::PathD p; for(auto&v:pts)p.emplace_back(v[0],v.y()); return p;
}
static std::vector<Vec2d> fromCP(const Clipper2Lib::PathD&p){
    std::vector<Vec2d> o; for(auto&v:p)o.emplace_back(v.x,v.y); return o;
}
Polygon2d SDFField::bufferPolygon(const Polygon2d&poly,double r){
    if(r<=0||poly.outer.empty())return poly;
    Clipper2Lib::PathsD ps; ps.push_back(toCP(poly.outer));
    auto inf=Clipper2Lib::InflatePaths(ps,r,Clipper2Lib::JoinType::Round,Clipper2Lib::EndType::Polygon);
    if(inf.empty())return poly;
    Polygon2d out; out.outer=fromCP(inf[0]); return out;
}
Polygon2d SDFField::unionPolygons(const std::vector<Polygon2d>&polys){
    if(polys.empty())return{};
    if(polys.size()==1)return polys[0];
    Clipper2Lib::PathsD sub;
    for(auto&p:polys)if(!p.outer.empty())sub.push_back(toCP(p.outer));
    auto sol=Clipper2Lib::Union(sub,Clipper2Lib::FillRule::NonZero);
    if(sol.empty())return polys[0];
    Polygon2d out; out.outer=fromCP(sol[0]); return out;
}
static double ptSegDist(const Vec2d&p,const Vec2d&a,const Vec2d&b){
    Vec2d ab=b-a,ap=p-a; double t=ab.dot(ap); double l2=ab.squaredNorm();
    if(l2<1e-20)return ap.norm(); t=std::max(0.0,std::min(1.0,t/l2));
    return(p-(a+t*ab)).norm();
}
double SDFField::distToPolygons(const Vec2d&pt,const std::vector<Polygon2d>&polys)const{
    double m=std::numeric_limits<double>::max();
    for(auto&poly:polys){auto&ring=poly.outer;int n=(int)ring.size();
        for(int i=0;i<n;++i)m=std::min(m,ptSegDist(pt,ring[i],ring[(i+1)%n]));}
    return m;
}
bool SDFField::insideAny(const Vec2d&pt,const std::vector<Polygon2d>&polys)const{
    for(auto&p:polys)if(polygonContains(p,pt))return true; return false;
}
void SDFField::build(const BoundingBox2d&roi,const std::vector<Obstacle>&obs,double cs,double buf){
    cs_=cs; roi_=roi;
    std::vector<Polygon2d> buffered;
    for(auto&o:obs){if(o.geometry.outer.empty())continue;
        buffered.push_back(buf>0?bufferPolygon(o.geometry,buf):o.geometry);}
    // FIX: do NOT force-union all buffered polygons into one.
    // When obstacles are non-overlapping, Clipper2's Union() returns multiple
    // separate paths and the old code only kept sol[0], silently discarding
    // all other obstacles from the SDF. Instead, keep every buffered polygon
    // independently; rebuildGrid / distToPolygons / insideAny already iterate
    // over all polygons in the list correctly.
    if(buffered.size()>1){
        // Attempt union only to merge OVERLAPPING obstacles (reduces polygon count).
        // For each connected component returned by Union, keep it as a separate
        // Polygon2d so non-overlapping obstacles are never discarded.
        Clipper2Lib::PathsD sub;
        for(auto&p:buffered)if(!p.outer.empty())sub.push_back(toCP(p.outer));
        auto sol=Clipper2Lib::Union(sub,Clipper2Lib::FillRule::NonZero);
        if(!sol.empty()){
            buffered_.clear();
            for(auto&path:sol){
                Polygon2d pg; pg.outer=fromCP(path); buffered_.push_back(pg);
            }
        } else {
            buffered_=buffered;
        }
    } else {
        buffered_=buffered;
    }
    rebuildGrid(buffered_);
}
void SDFField::buildFromPolygons(const BoundingBox2d&roi,const std::vector<Polygon2d>&polys,double cs){
    cs_=cs; roi_=roi; buffered_=polys; rebuildGrid(polys);
}
void SDFField::rebuildGrid(const std::vector<Polygon2d>&ops){
    double mg=cs_; BoundingBox2d ext;
    ext.min_pt=roi_.min_pt-Vec2d(mg,mg); ext.max_pt=roi_.max_pt+Vec2d(mg,mg); roi_=ext;
    cols_=std::max(2,(int)std::ceil(roi_.width()/cs_));
    rows_=std::max(2,(int)std::ceil(roi_.height()/cs_));
    grid_.assign(rows_*cols_,0.0);
    for(int r=0;r<rows_;++r)for(int c=0;c<cols_;++c){
        Vec2d pt=cellToWorld(r,c);
        double d=distToPolygons(pt,ops); bool ins=insideAny(pt,ops);
        grid_[idx(r,c)]=ins?-d:d;}
}
void SDFField::updateRegion(const BoundingBox2d&dirty,const std::vector<Obstacle>&obs,double buf){
    std::vector<Polygon2d> lp;
    for(auto&o:obs)lp.push_back(buf>0?bufferPolygon(o.geometry,buf):o.geometry);
    auto[rm,cm]=worldToCell(dirty.min_pt); auto[rx,cx]=worldToCell(dirty.max_pt);
    rm=std::max(0,rm-1);cm=std::max(0,cm-1);rx=std::min(rows_-1,rx+1);cx=std::min(cols_-1,cx+1);
    for(int r=rm;r<=rx;++r)for(int c=cm;c<=cx;++c){
        Vec2d pt=cellToWorld(r,c); double d=distToPolygons(pt,lp); bool ins=insideAny(pt,lp);
        grid_[idx(r,c)]=ins?-d:d;}
}
double SDFField::rawAt(int r,int c)const{
    r=std::max(0,std::min(rows_-1,r)); c=std::max(0,std::min(cols_-1,c)); return grid_[idx(r,c)];
}
std::pair<int,int> SDFField::worldToCell(const Vec2d&p)const{
    return{(int)((p[1]-roi_.min_pt.y())/cs_),(int)((p[0]-roi_.min_pt.x())/cs_)};
}
Vec2d SDFField::cellToWorld(int r,int c)const{
    return Vec2d(roi_.min_pt[0]+(c+0.5)*cs_,roi_.min_pt[1]+(r+0.5)*cs_);
}
std::pair<double,Vec2d> SDFField::queryWithGrad(const Vec2d&pt)const{
    if(grid_.empty())return{1e18,Vec2d(0,0)};
    double fx=(pt[0]-roi_.min_pt.x())/cs_-0.5,fy=(pt[1]-roi_.min_pt.y())/cs_-0.5;
    int x0=(int)std::floor(fx),y0=(int)std::floor(fy),x1=x0+1,y1=y0+1;
    double tx=fx-x0,ty=fy-y0;
    x0=std::max(0,std::min(cols_-1,x0));x1=std::max(0,std::min(cols_-1,x1));
    y0=std::max(0,std::min(rows_-1,y0));y1=std::max(0,std::min(rows_-1,y1));
    double v00=rawAt(y0,x0),v10=rawAt(y1,x0),v01=rawAt(y0,x1),v11=rawAt(y1,x1);
    double val=(1-tx)*(1-ty)*v00+tx*(1-ty)*v01+(1-tx)*ty*v10+tx*ty*v11;
    double gx=((1-ty)*(v01-v00)+ty*(v11-v10))/cs_;
    double gy=((1-tx)*(v10-v00)+tx*(v11-v01))/cs_;
    return{val,Vec2d(gx,gy)};
}
bool SDFField::isSafe(const Vec2d&pt,double cl)const{auto[d,_]=queryWithGrad(pt);return d>=cl;}
double SDFField::obstaclePenalty(const Vec2d&pt,double cl)const{
    auto[d,_]=queryWithGrad(pt);double s=d-cl;return s>=0?0.0:s*s;
}
Vec2d SDFField::obstaclePenaltyGrad(const Vec2d&pt,double cl)const{
    auto[d,gd]=queryWithGrad(pt);double s=d-cl;return s>=0?Vec2d(0,0):2.0*s*gd;
}
