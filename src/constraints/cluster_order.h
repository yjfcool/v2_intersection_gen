#pragma once

#include "types.h"
#include <unordered_map>

namespace isg {

enum class CrossExemption {
    None, StructuralCross, ObstacleCross
};

struct CurvePair {
    ConnId id_a, id_b;
    CrossExemption exempt = CrossExemption::None;
    double exempt_zone_radius = 0.0;
};

class SDFField;

class ClusterOrderSolver {
public:
    void build(const std::vector<Connectivity>&, const std::vector<Lane>&, const std::vector<LaneGroup>&);

    void markObstacleExempt(CurvePair&, const Vec2d&, const SDFField&, double r = 1.5);

    const std::vector<CurvePair>& pairs() const { return pairs_; }

    CrossExemption exemptionOf(const ConnId&, const ConnId&) const;

    void checkAndMarkA2(const std::unordered_map<ConnId, BezierCurve>&, const SDFField&, double r = 1.5);

    // Access to entry/exit lane orderings (for geometry-based tie-breaking in generator)
    const std::unordered_map<LaneId, std::vector<ConnId>>& entryOrder() const { return entry_order_; }
    const std::unordered_map<LaneId, std::vector<ConnId>>& exitOrder() const { return exit_order_; }

private:
    void addPairsFromCluster(const std::vector<ConnId>& cids, const std::vector<Connectivity>& conns);

    static bool hasPair(const std::vector<CurvePair>& pairs, const ConnId& a, const ConnId& b);

private:
    std::vector<CurvePair> pairs_;
    std::unordered_map<LaneId, std::vector<ConnId>> entry_order_;
    std::unordered_map<LaneId, std::vector<ConnId>> exit_order_;
};

}