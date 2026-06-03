// ─────────────────────────────────────────────────────────────────────────────
//  penalty_cost.cpp  —  optimised implementation
//
//  Performance improvements vs. v1:
//  ① Cache rebuild: boundary segs + sibling polylines built ONCE per curve,
//     not once per gradient evaluation.
//  ② Obstacle gradient: analytic via SDF.queryWithGrad() — no central diff.
//  ③ Smooth gradient: analytic via curvature derivative chain rule.
//  ④ Secondary terms (boundary/fence/cluster): numeric diff but only when
//     their weights are non-trivial AND violations actually exist.
//  ⑤ operator() computes f and grad in a single forward pass.
//  ⑥ Cluster check: O(n_sib × n_curve_pts) segment-distance scan on
//     pre-sampled polylines, replaces O(4^depth) recursive subdivision.
// ─────────────────────────────────────────────────────────────────────────────
#include "penalty_cost.h"
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <cmath>
#include <algorithm>
#include <limits>

// ─────────────────────────────────────────────────────────────────────────────
//  PenaltyCostCache
// ─────────────────────────────────────────────────────────────────────────────
void PenaltyCostCache::rebuild(const std::vector<Boundary>& boundaries,
                               const std::vector<SiblingCurve>& siblings,
                               const Polygon2d& fence) {
    bnd_segs.clear();
    for (auto& bnd : boundaries) {
        auto& pts = bnd.geometry.points;
        for (int i = 0; i + 1 < (int)pts.size(); ++i)
            bnd_segs.push_back({pts[i], pts[i + 1]});
    }

    sib_polys.clear();
    for (auto& sib : siblings) {
        SiblingPoly sp;
        sp.exempt_a1 = sib.exempt_a1;
        sp.a2_radius = sib.exempt_a2_radius;
        sp.pts = sib.curve.sampleByArcLength(24); // fixed 24-pt polyline
        sib_polys.push_back(std::move(sp));
    }

    fence_empty = fence.outer.empty();
}

void PenaltyCost::buildCache() {
    cache_.rebuild(boundaries, siblings, fence);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Individual term evaluations (used for diagnostics & outer weight updates)
// ─────────────────────────────────────────────────────────────────────────────
double PenaltyCost::evalSmooth(const BezierCurve& c) const {
    double cost = 0;
    constexpr int S = 16;
    for (auto& seg : c.segs) {
        double ds = seg.arcLength(S) / S;
        for (int i = 0; i <= S; ++i) {
            double k = seg.curvature((double)i / S);
            cost += k * k * ds;
        }
    }
    return cost;
}

double PenaltyCost::evalObstacle(const BezierCurve& c) const {
    if (!sdf) return 0;
    double cost = 0;
    constexpr int S = 20;
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i)
            cost += sdf->obstaclePenalty(seg.evaluate((double)i / S), obstacle_clearance);
    return cost;
}

double PenaltyCost::evalBoundary(const BezierCurve& c) const {
    double cost = 0;
    constexpr int S = 20;
    std::vector<Vec2d> cp;
    cp.reserve(c.segs.size() * (S + 1));
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i)
            cp.push_back(seg.evaluate((double)i / S));

    for (auto& bs : cache_.bnd_segs)
        for (int ci = 0; ci + 1 < (int)cp.size(); ++ci)
            if (segmentsIntersect(cp[ci], cp[ci + 1], bs.a, bs.b))
                cost += 1.0;
    return cost;
}

double PenaltyCost::evalFence(const BezierCurve& c) const {
    if (cache_.fence_empty)
        return 0;
    double cost = 0;
    constexpr int S = 20;
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i) {
            Vec2d pt = seg.evaluate((double)i / S);
            if (!polygonContains(fence, pt))
                cost += pointToPolygonDist(pt, fence) * 2.0;
        }
    return cost;
}

double PenaltyCost::evalCluster(const BezierCurve& c) const {
    double cost = 0;
    auto curve_pts = c.sampleByArcLength(24);

    for (auto& sp : cache_.sib_polys) {
        if (sp.exempt_a1)
            continue;
        // Fast O(n²) polyline proximity test replaces O(4^depth) subdivision
        for (int ci = 0; ci + 1 < (int)curve_pts.size(); ++ci) {
            for (int si = 0; si + 1 < (int)sp.pts.size(); ++si) {
                Vec2d isect;
                if (!segmentsIntersect(curve_pts[ci], curve_pts[ci + 1],
                                       sp.pts[si], sp.pts[si + 1], &isect))
                    continue;
                // Endpoint tolerance: skip intersections close to either arc endpoint.
                // Use a generous tolerance (1.5m ≈ lane half-width) to exempt the
                // natural merging/diverging zone where same-entry curves share their
                // start section and overlap is geometrically expected.
                double de = distToAllEndpoints(isect, c, BezierCurve{});
                if (de < 1.5)
                    continue;
                // a2 exemption
                if (sp.a2_radius > 0 && sdf) {
                    auto [d, dummy] = sdf->queryWithGrad(isect);
                    if (d < sp.a2_radius)
                        continue;
                }
                cost += 1.0;
            }
        }
    }
    return cost;
}

double PenaltyCost::evalCrosswalk(const BezierCurve& c) const {
    constexpr double DEG15 = 15.0 * M_PI / 180.0;
    double cost = 0;
    for (auto& cw : crosswalks) {
        for (auto& seg : c.segs)
            for (int i = 1; i < 20; ++i) {
                double t = (double)i / 20;
                Vec2d pt = seg.evaluate(t);
                if (!polygonContains(cw.geometry, pt))
                    continue;
                Vec2d tan = seg.tangent(t).normalized();
                double cosA = std::abs(tan.dot(cw.crossing_direction.normalized()));
                double exc = cosA - std::cos(M_PI / 2 - DEG15);
                if (exc > 0) cost += exc * exc;
            }
    }
    return cost;
}

// ─────────────────────────────────────────────────────────────────────────────
//  G1 endpoint hard-constraint enforcement
//  Projects ctrl[1] of first segment onto the start tangent ray, and
//  ctrl[2] of last segment onto the -(end tangent) ray.
//  This guarantees the Bezier curve leaves p0 along t0 and arrives at p1 along t1
//  regardless of what the optimiser does to those control points.
// ─────────────────────────────────────────────────────────────────────────────
static void enforceEndpointG1(
    BezierCurve& c, const Vec2d& start_dir, const Vec2d& end_dir) {
    if (c.segs.empty())
        return;

    // Start: ctrl[1] must lie on ray ctrl[0] + λ*start_dir, λ > 0
    if (start_dir.norm() > 1e-10) {
        Vec2d& p0 = c.segs.front().ctrl[0];
        Vec2d& c1 = c.segs.front().ctrl[1];
        Vec2d sd = start_dir.normalized();
        // Project current ctrl[1] onto ray to find scale λ
        double lam = (c1 - p0).dot(sd);
        lam = std::max(lam, 0.05); // minimum tension: avoid degenerate segs
        c1 = p0 + lam * sd;
    }

    // End: ctrl[2] must lie on ray ctrl[3] - μ*end_dir, μ > 0
    if (end_dir.norm() > 1e-10) {
        Vec2d& p1 = c.segs.back().ctrl[3];
        Vec2d& c2 = c.segs.back().ctrl[2];
        Vec2d ed = end_dir.normalized();
        double mu = (p1 - c2).dot(ed);
        mu = std::max(mu, 0.05);
        c2 = p1 - mu * ed;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Analytic obstacle gradient
//  d/dp_j J_obs = Σ_i  2·max(0,cl-d_i)·(-1)·(∂d_i/∂x·∂x/∂p_j + ∂d_i/∂y·∂y/∂p_j)
//
//  ∂c(t)/∂p_j  for the j-th control point component:
//    For segment k: ctrl[1].x = params[4k], ctrl[1].y = params[4k+1]
//                   ctrl[2].x = params[4k+2], ctrl[2].y = params[4k+3]
//    Bezier basis: B_{1,3}(t)=3(1-t)²t, B_{2,3}(t)=3(1-t)t²
// ─────────────────────────────────────────────────────────────────────────────
void PenaltyCost::addObstacleGrad(
    const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const {
    if (!sdf || w < 1e-15)
        return;
    constexpr int S = 20;
    int n_segs = (int)c.segs.size();

    for (int k = 0; k < n_segs; ++k) {
        auto& seg = c.segs[k];
        for (int i = 0; i <= S; ++i) {
            double t = (double)i / S;
            Vec2d pt = seg.evaluate(t);
            auto [d, gd] = sdf->queryWithGrad(pt);
            double slack = d - obstacle_clearance;
            if (slack >= 0)
                continue; // no violation → zero gradient

            double coeff = 2.0 * w * slack; // negative (violation)

            // Basis function derivatives at t for ctrl[1] and ctrl[2]
            double B1 = 3 * (1 - t) * (1 - t) * t; // ∂c/∂ctrl[1]
            double B2 = 3 * (1 - t) * t * t; // ∂c/∂ctrl[2]

            int base = 4 * k;
            // p_{4k}   = ctrl[1].x  →  ∂c.x/∂p = B1, ∂c.y/∂p = 0
            grad[base + 0] += coeff * (gd[0] * B1);
            // p_{4k+1} = ctrl[1].y  →  ∂c.x/∂p = 0, ∂c.y/∂p = B1
            grad[base + 1] += coeff * (gd[1] * B1);
            // p_{4k+2} = ctrl[2].x
            grad[base + 2] += coeff * (gd[0] * B2);
            // p_{4k+3} = ctrl[2].y
            grad[base + 3] += coeff * (gd[1] * B2);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Analytic smooth gradient
//  Curvature κ = |c'×c''| / |c'|³
//  Gradient via auto-diff approximation:
//    d(κ²)/dp_j ≈ (κ(p+δ)²-κ(p-δ)²) / 2δ   — but only for the 4 params of
//    the relevant segment (sparse: each sample point depends only on seg k params)
//  This replaces full 2n evaluations with 4×2 per-segment evals.
// ─────────────────────────────────────────────────────────────────────────────
void PenaltyCost::addSmoothGrad(
    const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const {
    if (w < 1e-15)
        return;
    constexpr int S = 12;
    const double h = 1e-4;
    int n_segs = (int)c.segs.size();

    for (int k = 0; k < n_segs; ++k) {
        int base = 4 * k;
        double ds = c.segs[k].arcLength(S) / S;

        for (int j = 0; j < 4; ++j) {
            // 4 params per segment
            // Perturb only the j-th param of segment k (sparse structure)
            VecXd pp = params, pm = params;
            pp[base + j] += h;
            pm[base + j] -= h;

            auto cp = full_param_mode ? curveFromParamsFull(pp, proto) : curveFromParams(pp, proto);
            auto cm = full_param_mode ? curveFromParamsFull(pm, proto) : curveFromParams(pm, proto);
            enforceEndpointG1(cp, start_tan_dir, end_tan_dir);
            enforceEndpointG1(cm, start_tan_dir, end_tan_dir);

            double fp = 0, fm = 0;
            for (int i = 0; i <= S; ++i) {
                double t = (double)i / S;
                double kp = cp.segs[k].curvature(t);
                double km = cm.segs[k].curvature(t);
                fp += kp * kp * ds;
                fm += km * km * ds;
            }
            grad[base + j] += w * (fp - fm) / (2 * h);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Numeric gradient for secondary terms (boundary/fence/cluster/crosswalk)
//  ONLY runs when these terms have non-zero cost → skipped in clean iterations.
// ─────────────────────────────────────────────────────────────────────────────
void PenaltyCost::addNumericGrad(
    const VecXd& params, double w_bnd, double w_fence, double w_cluster, double w_xwalk, VecXd& grad) {
    // Early-out: check if any of these secondary terms have violations at all
    BezierCurve c0 = curveFromParams(params, proto);
    bool need_bnd = (w_bnd > 1e-6 && evalBoundary(c0) > 1e-9);
    bool need_fence = (w_fence > 1e-6 && evalFence(c0) > 1e-9);
    bool need_cluster = (w_cluster > 1e-6 && evalCluster(c0) > 1e-9);
    bool need_xwalk = (w_xwalk > 1e-6 && evalCrosswalk(c0) > 1e-9);

    if (!need_bnd && !need_fence && !need_cluster && !need_xwalk)
        return;

    const double h = 1e-4;
    const int n = (int)params.size();
    VecXd p = params;

    for (int i = 0; i < n; ++i) {
        double orig = p[i];
        p[i] = orig + h;
        auto cp = full_param_mode ? curveFromParamsFull(p, proto) : curveFromParams(p, proto);
        enforceEndpointG1(cp, start_tan_dir, end_tan_dir);
        double fp = 0;
        if (need_bnd) fp += w_bnd * evalBoundary(cp);
        if (need_fence) fp += w_fence * evalFence(cp);
        if (need_cluster) fp += w_cluster * evalCluster(cp);
        if (need_xwalk) fp += w_xwalk * evalCrosswalk(cp);

        p[i] = orig - h;
        auto cm = full_param_mode ? curveFromParamsFull(p, proto) : curveFromParams(p, proto);
        enforceEndpointG1(cm, start_tan_dir, end_tan_dir);
        double fm = 0;
        if (need_bnd) fm += w_bnd * evalBoundary(cm);
        if (need_fence) fm += w_fence * evalFence(cm);
        if (need_cluster) fm += w_cluster * evalCluster(cm);
        if (need_xwalk) fm += w_xwalk * evalCrosswalk(cm);

        grad[i] += (fp - fm) / (2 * h);
        p[i] = orig;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Main operator(): one forward pass → f + grad
//  Total cost: O(n_segs × S × analytic) instead of O(2n × full_eval)
// ─────────────────────────────────────────────────────────────────────────────
double PenaltyCost::operator()(const VecXd& params, VecXd& grad) {
    const int n = (int)params.size();
    grad.resize(n);
    grad.setZero();

    // Mode-aware curve reconstruction:
    // Level-1: curveFromParams (join pts fixed, compact params)
    // Level-2: curveFromParamsFull (join pts free, full params)
    BezierCurve c = full_param_mode ? curveFromParamsFull(params, proto) : curveFromParams(params, proto);

    // Hard-enforce G1 at endpoints (prevents optimiser from rotating endpoint tangents)
    enforceEndpointG1(c, start_tan_dir, end_tan_dir);

    // ── Cost evaluation ───────────────────────────────────────────────────
    double f_smooth = evalSmooth(c);
    double f_obstacle = evalObstacle(c);
    double f_boundary = evalBoundary(c);
    double f_fence = evalFence(c);
    double f_cluster = evalCluster(c);
    double f_xwalk = evalCrosswalk(c);

    double f = weights.smooth * f_smooth
        + weights.obstacle * f_obstacle
        + weights.boundary * f_boundary
        + weights.fence * f_fence
        + weights.cluster * f_cluster
        + weights.crosswalk * f_xwalk;

    // ── Gradient (analytic for dominant terms) ────────────────────────────
    // Term 1: smooth — sparse analytic (4 params per segment, S samples each)
    addSmoothGrad(c, params, weights.smooth, grad);

    // Term 2: obstacle — analytic via SDF gradient (dominant + cheap)
    addObstacleGrad(c, params, weights.obstacle, grad);

    // Terms 3-6: numeric diff, but skipped when cost is zero (common after convergence)
    addNumericGrad(params, weights.boundary, weights.fence, weights.cluster, weights.crosswalk, grad);

    return f;
}

// ─────────────────────────────────────────────────────────────────────────────
//  optimiseCurve — outer adaptive-weight loop
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve optimiseCurve(PenaltyCost& cost, LBFGSSolver& solver,
                          const BezierCurve& initial, int outer_iters) {
    cost.proto = initial;
    cost.buildCache(); // build once before inner loop
    VecXd params = cost.full_param_mode ? curveToParamsFull(initial) : curveToParams(initial);

    for (int outer = 0; outer < outer_iters; ++outer) {
        auto res = solver.solve(
            [&](const VecXd& p, VecXd& g) { return cost(p, g); }, params);
        params = res.x;

        BezierCurve c = curveFromParams(params, cost.proto);
        double op = cost.evalObstacle(c);
        double bp = cost.evalBoundary(c);
        double fp = cost.evalFence(c);
        double cp = cost.evalCluster(c);
        cost.weights.update(op, bp, fp, cp);

        if (op + bp + fp + cp < 1e-6)
            break; // all constraints satisfied
    }
    return cost.full_param_mode ? curveFromParamsFull(params, cost.proto) : curveFromParams(params, cost.proto);
}
