#include "cluster_order.h"
#include "optimizer/sdf_field.h"
#include <algorithm>

#include "curve/curve_utils.h"

// ─────────────────────────────────────────────────────────────────────────────
//  Turn-type geometry order for non-crossing constraint within a cluster.
//
//  Correct inner→outer ordering based on actual arc geometry at a junction:
//    Left U-turn  (most inner on left side)
//    Left turn
//    Straight     (reference / baseline)
//    Right turn
//    Right U-turn (most inner on right side, but outer arc for the cluster)
//
//  The key insight: "inner" means smaller arc radius / shorter path.
//  - Left-turn curves arc further left (inner for the left side cluster),
//    so they should be sorted innermost when leaving from the same entry.
//  - Straight is the middle baseline.
//  - Right-turn curves arc right (inner for the right side cluster).
//  - U-turns are large semicircular arcs that extend far, so they occupy
//    the outermost position (they swing furthest from the straight baseline).
//
//  The ordering priority below encodes "distance from straight baseline",
//  with straight=0 (reference), turns stepping away, U-turns furthest out.
//  Within same priority, multiple connections to different exit lanes are
//  sorted geometrically (see sortSameTypeBySpatial).
//
//  IMPORTANT: UTurn exemption REMOVED per requirements. All pairs within
//  a cluster (same entry or same exit) must respect non-crossing, unless
//  there is a genuine logical conflict (different road arms crossing).
// ─────────────────────────────────────────────────────────────────────────────

// Priority assigns a "lateral offset category" relative to straight:
//   UTurnLeft  → -2 (swings hard left, innermost from left perspective)
//   TurnLeft   → -1
//   Straight   →  0 (baseline)
//   TurnRight  → +1
//   UTurnRight → +2
// Sort order: by ascending priority (left first → right last).
static int turnPriority(ConnTurnType t) {
    switch (t) {
    case ConnTurnType::UTurnLeft:  return -2;
    case ConnTurnType::TurnLeft:   return -1;
    case ConnTurnType::Straight:   return  0;
    case ConnTurnType::TurnRight:  return +1;
    case ConnTurnType::UTurnRight: return +2;
    default:                        return  0;
    }
}

// Returns the Connectivity for a given id (nullptr if not found).
static const Connectivity* findConn(const std::vector<Connectivity>& conns, const ConnId& id) {
    for (auto& c : conns)
        if (c.id == id) return &c;
    return nullptr;
}

// Check if a pair (a, b) or (b, a) already exists in the list.
bool ClusterOrderSolver::hasPair(const std::vector<CurvePair>& pairs,
                                 const ConnId& a, const ConnId& b) {
    for (auto& p : pairs)
        if ((p.id_a == a && p.id_b == b) || (p.id_a == b && p.id_b == a))
            return true;
    return false;
}

// Add non-crossing pairs for a cluster (entry or exit).
// Pairs within the same cluster should be ordered by turn priority.
// No automatic UTurn exemption — let the optimizer handle it or mark
// topological conflicts via checkAndMarkA2 / infeasibility detection.
void ClusterOrderSolver::addPairsFromCluster(
    const std::vector<ConnId>& cids,
    const std::vector<Connectivity>& conns)
{
    for (int i = 0; i < (int)cids.size(); ++i) {
        for (int j = i + 1; j < (int)cids.size(); ++j) {
            if (hasPair(pairs_, cids[i], cids[j]))
                continue;
            CurvePair p;
            p.id_a = cids[i];
            p.id_b = cids[j];
            // No automatic exemption — all same-cluster pairs are constrained.
            // (Topological conflicts can be detected later by checkAndMarkA2.)
            pairs_.push_back(p);
        }
    }
}

void ClusterOrderSolver::build(const std::vector<Connectivity>& conns,
                               const std::vector<Lane>& lanes,
                               const std::vector<LaneGroup>& laneGroups) {
    entry_order_.clear();
    exit_order_.clear();
    pairs_.clear();

    // ── Step 1: Group connections by entry lane and exit lane ────────────────
    for (auto& c : conns) {
        entry_order_[c.entry_lane_id].push_back(c.id);
        exit_order_ [c.exit_lane_id ].push_back(c.id);
    }

    // ── Step 2: Sort each entry cluster by turn-type priority ────────────────
    // Sort rule: ascending turnPriority (UTurnLeft < TurnLeft < Straight < TurnRight < UTurnRight)
    // Within same priority, keep insertion order (stable, geometry-preserving).
    // For same-turn-type multiple connections (e.g. 3 right turns to 3 different exits),
    // we sort by exit lane identity to get a deterministic consistent ordering.
    // The generator will further refine same-type ordering geometrically.
    for (auto& kv : entry_order_) {
        auto& lid = kv.first; (void)lid; auto& cids = kv.second;
        std::stable_sort(cids.begin(), cids.end(), [&](const ConnId& a, const ConnId& b) {
            auto* ca = findConn(conns, a);
            auto* cb = findConn(conns, b);
            ConnTurnType ta = ca ? ca->turn_type : ConnTurnType::Straight;
            ConnTurnType tb = cb ? cb->turn_type : ConnTurnType::Straight;
            int pa = turnPriority(ta), pb = turnPriority(tb);
            if (pa != pb) return pa < pb;
            // Same turn type: sort by exit_lane_id for stable ordering
            LaneId xa = ca ? ca->exit_lane_id : "";
            LaneId xb = cb ? cb->exit_lane_id : "";
            return xa < xb;
        });
    }

    // ── Step 3: Sort each exit cluster by turn-type priority (reversed) ──────
    // For a same-exit cluster, the natural order is OPPOSITE to entry order:
    // curves arriving from different directions at the same exit.
    // Use the same priority scheme — let the optimizer handle ordering.
    for (auto& kv : exit_order_) {
        auto& lid = kv.first; (void)lid; auto& cids = kv.second;
        std::stable_sort(cids.begin(), cids.end(), [&](const ConnId& a, const ConnId& b) {
            auto* ca = findConn(conns, a);
            auto* cb = findConn(conns, b);
            ConnTurnType ta = ca ? ca->turn_type : ConnTurnType::Straight;
            ConnTurnType tb = cb ? cb->turn_type : ConnTurnType::Straight;
            int pa = turnPriority(ta), pb = turnPriority(tb);
            if (pa != pb) return pa < pb;
            LaneId ea = ca ? ca->entry_lane_id : "";
            LaneId eb = cb ? cb->entry_lane_id : "";
            return ea < eb;
        });
    }

    // ── Step 4: Build constraint pairs from entry clusters ───────────────────
    for (auto& kv : entry_order_) {
        auto& lid = kv.first; (void)lid; auto& cids = kv.second;
        if (cids.size() > 1)
            addPairsFromCluster(cids, conns);
    }

    // ── Step 5: Build constraint pairs from exit clusters ────────────────────
    // This is the KEY FIX: curves sharing the same exit lane must not cross
    // (e.g. conn32/right-turn and conn18/left-turn sharing exit 43242048).
    for (auto& kv : exit_order_) {
        auto& lid = kv.first; (void)lid; auto& cids = kv.second;
        if (cids.size() > 1)
            addPairsFromCluster(cids, conns);
    }
}

void ClusterOrderSolver::markObstacleExempt(CurvePair& pair, const Vec2d& pt, const SDFField& sdf, double r) {
    std::pair<double,Vec2d> _q = sdf.queryWithGrad(pt);
    if (_q.first < r) {
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
