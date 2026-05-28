#pragma once
#include "types.h"
#include "sdf_field.h"
#include "lbfgs_solver.h"
#include <vector>

struct SiblingCurve {
    BezierCurve curve;
    std::vector<Vec2d> sample_pts; // pre-sampled for fast cluster check
    bool exempt_a1 = false;
    double exempt_a2_radius = 0.0;
};

struct PenaltyWeights {
    double smooth = 1.0;
    double obstacle = 10.0;
    double boundary = 8.0;
    double fence = 5.0;
    double cluster = 3.0;
    double crosswalk = 0.5;

    void update(double op, double bp, double fp, double cp) {
        if (op > 1e-3) obstacle = std::min(obstacle * 2.0, 160.0);
        if (bp > 1e-3) boundary = std::min(boundary * 2.0, 128.0);
        if (fp > 1e-3) fence = std::min(fence * 2.0, 80.0);
        if (cp > 1e-3) cluster = std::min(cluster * 1.5, 48.0);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  Pre-computed spatial data (built once per PenaltyCost instance, not per eval)
// ─────────────────────────────────────────────────────────────────────────────
struct BoundarySegment {
    Vec2d a, b;
};

struct PenaltyCostCache {
    // Per-sibling pre-sampled polylines (avoids subdivision tree on every call)
    // Updated when siblings change; accessed read-only during gradient computation.
    struct SiblingPoly {
        std::vector<Vec2d> pts;
        bool exempt_a1;
        double a2_radius;
    };

    // All boundary line segments in one flat list (avoid re-iterating per eval)
    std::vector<BoundarySegment> bnd_segs;
    std::vector<SiblingPoly> sib_polys;
    // Fence polygon winding-number lookup cached at construction
    bool fence_empty = true;

    void rebuild(const std::vector<Boundary>& bounds, const std::vector<SiblingCurve>& sibls, const Polygon2d& fence);
};

// ─────────────────────────────────────────────────────────────────────────────
//  PenaltyCost — unified objective function
//
//  Performance highlights vs. v1:
//  1. Analytic gradient for obstacle & smooth terms → 16× fewer full evals
//  2. Pre-sampled sibling polylines → cluster check is O(n) not O(4^depth)
//  3. Pre-flattened boundary segment list → boundary check is O(n_bnd × n_curve)
//  4. curveFromParams operates in-place on a scratch BezierCurve
//  5. operator() computes f + analytic grad in ONE pass, no central differences
// ─────────────────────────────────────────────────────────────────────────────
class PenaltyCost {
public:
    BezierCurve proto;
    PenaltyWeights weights;
    const SDFField* sdf = nullptr;
    std::vector<Boundary> boundaries;
    Polygon2d fence;
    std::vector<SiblingCurve> siblings;
    std::vector<Crosswalk> crosswalks;
    double obstacle_clearance = 0.1;

    // Level-2 mode: use curveToParamsFull/curveFromParamsFull
    // so join points are also optimised (no frozen topology)
    bool full_param_mode = false;

    // G1 hard-constraint: endpoint tangent directions (unit vectors)
    // If set (non-zero), curveFromParams() enforces ctrl[1] along start_tan
    // and ctrl[2] of last seg along -end_tan, preserving G1 at all times.
    Vec2d start_tan_dir{0, 0}; // set to t0 from Lane::geometry
    Vec2d end_tan_dir{0, 0}; // set to t1 from Lane::geometry

    // Must be called once after setting boundaries/siblings/fence
    void buildCache();

    // Main entry point: computes f and analytic+numeric-hybrid gradient
    double operator()(const VecXd& params, VecXd& grad);

    // Individual term evaluations (for diagnostics & weight update)
    double evalSmooth(const BezierCurve& c) const;
    double evalObstacle(const BezierCurve& c) const;
    double evalBoundary(const BezierCurve& c) const;
    double evalFence(const BezierCurve& c) const;
    double evalCluster(const BezierCurve& c) const;
    double evalCrosswalk(const BezierCurve& c) const;

private:
    PenaltyCostCache cache_;

    // Compute analytic gradient contributions from obstacle SDF term
    // d/dp_j [ Σ_i max(0, cl - d(c(t_i)))² ]
    void addObstacleGrad(const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const;

    // Numeric central-diff gradient for boundary+fence+cluster (cheap compared
    // to obstacle since these terms are usually zero after early iterations)
    void addNumericGrad(const VecXd& params, double w_bnd, double w_fence, double w_cluster, double w_xwalk,
                        VecXd& grad);

    // Smooth term analytic gradient via curvature derivative
    void addSmoothGrad(const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const;
};

BezierCurve optimiseCurve(PenaltyCost& cost, LBFGSSolver& solver, const BezierCurve& initial, int outer_iters = 4);
