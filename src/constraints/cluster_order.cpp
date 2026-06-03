#include "cluster_order.h"
#include "optimizer/sdf_field.h"
#include <algorithm>

#include "curve/curve_utils.h"

void ClusterOrderSolver::build(const std::vector<Connectivity>& conns,
                               const std::vector<Lane>&, const std::vector<LaneGroup>&) {
    // Turn type priority for inner→outer ordering within a cluster.
    // (UTurnLeft=0 is most inner, UTurnRight=4 most outer)
    auto pri = [](ConnTurnType t) {
        switch (t) {
        case ConnTurnType::UTurnLeft:  return 0;
        case ConnTurnType::TurnLeft:   return 1;
        case ConnTurnType::Straight:   return 2;
        case ConnTurnType::TurnRight:  return 3;
        case ConnTurnType::UTurnRight: return 4;
        }
        return 2;
    };
    for (auto& c : conns)
        entry_order_[c.entry_lane_id].push_back(c.id);
    for (auto& [lid,cids] : entry_order_) {
        std::sort(cids.begin(), cids.end(), [&](const ConnId& a, const ConnId& b) {
            ConnTurnType ta = ConnTurnType::Straight, tb = ConnTurnType::Straight;
            for (auto& c : conns) {
                if (c.id == a)
                    ta = c.turn_type;
                if (c.id == b)
                    tb = c.turn_type;
            }
            return pri(ta) < pri(tb);
        });
    }
    pairs_.clear();
    for (auto& [lid,cids] : entry_order_) {
        for (int i = 0; i < (int)cids.size(); ++i)
            for (int j = i + 1; j < (int)cids.size(); ++j) {
                CurvePair p;
                p.id_a = cids[i];
                p.id_b = cids[j];
                ConnTurnType ta = ConnTurnType::Straight, tb = ConnTurnType::Straight;
                for (auto& c : conns) {
                    if (c.id == p.id_a) ta = c.turn_type;
                    if (c.id == p.id_b) tb = c.turn_type;
                }
                bool au = (ta == ConnTurnType::UTurnLeft || ta == ConnTurnType::UTurnRight);
                bool bu = (tb == ConnTurnType::UTurnLeft || tb == ConnTurnType::UTurnRight);
                if (au || bu)
                    p.exempt = CrossExemption::StructuralCross;
                pairs_.push_back(p);
            }
    }
}

void ClusterOrderSolver::markObstacleExempt(CurvePair& pair, const Vec2d& pt, const SDFField& sdf, double r) {
    auto [d,_] = sdf.queryWithGrad(pt);
    if (d < r) {
        pair.exempt = CrossExemption::ObstacleCross;
        pair.exempt_zone_radius = r;
    }
}

CrossExemption ClusterOrderSolver::exemptionOf(const ConnId& a, const ConnId& b) const {
    for (auto& p : pairs_)
        if ((p.id_a == a && p.id_b == b) || (p.id_a == b && p.id_b == a))
            return p.exempt;
    return CrossExemption::None;
}

void ClusterOrderSolver::checkAndMarkA2(
    const std::unordered_map<ConnId, BezierCurve>& curves, const SDFField& sdf, double r) {
    for (auto& pair : pairs_) {
        if (pair.exempt != CrossExemption::None)
            continue;
        auto ia = curves.find(pair.id_a), ib = curves.find(pair.id_b);
        if (ia == curves.end() || ib == curves.end())
            continue;
        for (auto& pt : curveCrossings(ia->second, ib->second))
            markObstacleExempt(pair, pt, sdf, r);
    }
}
