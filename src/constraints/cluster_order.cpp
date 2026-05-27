#include "cluster_order.h"
#include "optimizer/sdf_field.h"
#include "constraints/intersection_check.h"
#include <algorithm>
void ClusterOrderSolver::build(const std::vector<Connectivity>&conns,
    const std::vector<Lane>&,const std::vector<LaneGroup>&)
{
    auto pri=[](TurnType t){
        switch(t){case TurnType::UTurnLeft:return 0;case TurnType::TurnLeft:return 1;
        case TurnType::Straight:return 2;case TurnType::TurnRight:return 3;case TurnType::UTurnRight:return 4;}return 2;};
    for(auto&c:conns)entry_order_[c.entry_lane_id].push_back(c.id);
    for(auto&[lid,cids]:entry_order_){
        std::sort(cids.begin(),cids.end(),[&](const ConnId&a,const ConnId&b){
            TurnType ta=TurnType::Straight,tb=TurnType::Straight;
            for(auto&c:conns){if(c.id==a)ta=c.turn_type;if(c.id==b)tb=c.turn_type;}
            return pri(ta)<pri(tb);});
    }
    pairs_.clear();
    for(auto&[lid,cids]:entry_order_){
        for(int i=0;i<(int)cids.size();++i)for(int j=i+1;j<(int)cids.size();++j){
            CurvePair p; p.id_a=cids[i]; p.id_b=cids[j];
            TurnType ta=TurnType::Straight,tb=TurnType::Straight;
            for(auto&c:conns){if(c.id==p.id_a)ta=c.turn_type;if(c.id==p.id_b)tb=c.turn_type;}
            bool au=(ta==TurnType::UTurnLeft||ta==TurnType::UTurnRight);
            bool bu=(tb==TurnType::UTurnLeft||tb==TurnType::UTurnRight);
            if(au||bu)p.exempt=CrossExemption::StructuralCross;
            pairs_.push_back(p);
        }
    }
}
void ClusterOrderSolver::markObstacleExempt(CurvePair&pair,const Vec2d&pt,const SDFField&sdf,double r){
    auto[d,_]=sdf.queryWithGrad(pt);
    if(d<r){pair.exempt=CrossExemption::ObstacleCross;pair.exempt_zone_radius=r;}
}
CrossExemption ClusterOrderSolver::exemptionOf(const ConnId&a,const ConnId&b)const{
    for(auto&p:pairs_)if((p.id_a==a&&p.id_b==b)||(p.id_a==b&&p.id_b==a))return p.exempt;
    return CrossExemption::None;
}
void ClusterOrderSolver::checkAndMarkA2(
    const std::unordered_map<ConnId,BezierCurve>&curves,const SDFField&sdf,double r)
{
    for(auto&pair:pairs_){
        if(pair.exempt!=CrossExemption::None)continue;
        auto ia=curves.find(pair.id_a),ib=curves.find(pair.id_b);
        if(ia==curves.end()||ib==curves.end())continue;
        for(auto&pt:curveCrossings(ia->second,ib->second))
            markObstacleExempt(pair,pt,sdf,r);
    }
}
