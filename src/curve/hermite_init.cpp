// ─────────────────────────────────────────────────────────────────────────────
//  hermite_init.cpp  —  Two-level obstacle avoidance strategy
//
//  Level-1 (方案A): Geometric Direct Construction
//    Fast-path for ≤3 convex obstacles.
//    Steps: straightLine check → obstacle AABB union → compute apex
//    → build 2-segment G1 Bezier arch.  Zero iterations, analytic quality.
//
//  Level-2 (方案C): Full control-point optimisation
//    Fallback for multi-obstacle / concave / Level-1 failure.
//    ALL control points (incl. join points) are optimisation variables.
//    No frozen topology. SDF penalty drives smooth arch naturally.
// ─────────────────────────────────────────────────────────────────────────────
#include "hermite_init.h"
#include "optimizer/sdf_field.h"
#include "constraints/fence_check.h"
#include <cmath>
#include <algorithm>
#include <limits>

// ─────────────────────────────────────────────────────────────────────────────
//  Shared helpers
// ─────────────────────────────────────────────────────────────────────────────
static bool segmentsIntersect_internal(const Vec2d& a, const Vec2d& b,
                                        const Vec2d& c, const Vec2d& d,
                                        Vec2d* out = nullptr)
{
    Vec2d r = b-a, s = d-c;
    double den = cross2d(r, s);
    if (std::abs(den) < 1e-12) return false;
    Vec2d ac = c-a;
    double t = cross2d(ac, s) / den;
    double u = cross2d(ac, r) / den;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        if (out) *out = a + t*r;
        return true;
    }
    return false;
}
static constexpr double MAX_TAN_DEV = 60.0 * M_PI / 180.0;

static Vec2d clampTangent(const Vec2d& tan, const Vec2d& ref, double max_a) {
    Vec2d t = tan.norm() > 1e-10 ? tan.normalized() : ref;
    if (t.dot(ref) >= std::cos(max_a)) return t;
    double sg = cross2d(ref, t) >= 0 ? 1.0 : -1.0;
    double ca = std::cos(max_a), sa = std::sin(max_a);
    return Vec2d(ref.x()*ca - sg*ref.y()*sa,
                 ref.x()*sa*sg + ref.y()*ca);
}

// SDF-sampled straight line clearance check
static bool straightLineClear(const SDFField& sdf,
                               const Vec2d& p0, const Vec2d& p1,
                               double clearance = 0.15, int n = 30)
{
    if (!sdf.valid()) return true;
    for (int i = 1; i < n; ++i) {
        double t = (double)i / n;
        auto [d, dummy] = sdf.queryWithGrad((1-t)*p0 + t*p1);
        if (d < clearance) return false;
    }
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Level-1: Geometric Direct Construction (方案A)
// ─────────────────────────────────────────────────────────────────────────────

// Union-AABB of all obstacle buffered geometries visible to the SDF field.
// We probe the SDF along the direct path to find the obstacle's bounding box
// without needing obstacle polygon data (SDF is already available).
struct ObstacleAABB {
    double x_min, x_max, y_min, y_max;
    bool valid = false;
};

static ObstacleAABB probeObstacleAABB(const SDFField& sdf,
                                        const Vec2d& p0, const Vec2d& p1,
                                        double clearance = 0.15)
{
    // Sample dense grid in the neighbourhood of the direct path
    ObstacleAABB box;
    box.x_min = box.y_min =  1e18;
    box.x_max = box.y_max = -1e18;

    int nx = 60, ny = 30;
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) return box;
    along = along * (1.0/len);
    Vec2d perp{-along.y(), along.x()};
    double sweep = std::max(3.0, 0.4 * len);   // lateral scan width

    for (int ix = 0; ix <= nx; ++ix) {
        for (int iy = -ny/2; iy <= ny/2; ++iy) {
            double s = (double)ix / nx * len;
            double t = (double)iy / (ny/2) * sweep;
            Vec2d pt = p0 + s*along + t*perp;
            auto [d, dummy] = sdf.queryWithGrad(pt);
            if (d < clearance) {   // in or near obstacle
                box.x_min = std::min(box.x_min, pt.x());
                box.x_max = std::max(box.x_max, pt.x());
                box.y_min = std::min(box.y_min, pt.y());
                box.y_max = std::max(box.y_max, pt.y());
                box.valid = true;
            }
        }
    }
    return box;
}

// Compute the bypass apex point for a given side (+1 = above, -1 = below).
// Apex is placed at the obstacle AABB edge + safe clearance, at the
// longitudinal midpoint of the obstacle.
static Vec2d computeApex(const ObstacleAABB& box,
                           const Vec2d& p0, const Vec2d& p1,
                           int side,              // +1 = left/above, -1 = right/below
                           double clearance)
{
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) return 0.5*(p0+p1);
    along = along * (1.0/len);
    Vec2d perp{-along.y(), along.x()};   // left normal

    // Longitudinal midpoint of obstacle
    double lon_mid = 0.5*(box.x_min + box.x_max) * along.x()
                   + 0.5*(box.y_min + box.y_max) * along.y();
    // Clamp to [20%, 80%] of path length
    double path_lon_0 = p0.dot(along);
    double path_lon_1 = p1.dot(along);
    double frac = (lon_mid - path_lon_0) / (path_lon_1 - path_lon_0 + 1e-12);
    frac = std::max(0.2, std::min(0.8, frac));

    // Lateral position of apex: obstacle edge + clearance, on chosen side
    double lat_edge;
    if (side > 0)
        lat_edge = box.y_max * perp.y() + box.x_max * perp.x() + clearance;
    else
        lat_edge = box.y_min * perp.y() + box.x_min * perp.x() - clearance;

    // Apex in world coordinates
    Vec2d base = p0 + frac * (p1 - p0);   // point on direct line at lon_frac
    // Project lateral: apex.dot(perp) = lat_edge
    double base_lat = base.dot(perp);
    Vec2d apex = base + (lat_edge - base_lat) * perp;
    return apex;
}

// Choose which side (left/right) gives more clearance.
// Also receives already-generated sibling curves to avoid picking
// the same side as an opposing-direction arch (prevents intersection).
static int chooseSide(const SDFField& sdf,
                       const ObstacleAABB& box,
                       const Vec2d& p0, const Vec2d& p1,
                       double clearance,
                       const std::vector<std::vector<Vec2d>>& sibling_polys)
{
    Vec2d along = (p1 - p0).normalized();
    Vec2d perp{-along.y(), along.x()};

    double lat_spread = std::max(1.5,
        std::max(std::abs(box.y_max), std::abs(box.y_min)) + clearance + 0.5);

    Vec2d test_above = 0.5*(p0+p1) + lat_spread * perp;
    Vec2d test_below = 0.5*(p0+p1) - lat_spread * perp;

    auto [d_above, dum1] = sdf.queryWithGrad(test_above);
    auto [d_below, dum2] = sdf.queryWithGrad(test_below);

    // Count how many sibling sample points occupy each side
    int sib_above = 0, sib_below = 0;
    for (auto& poly : sibling_polys) {
        for (auto& pt : poly) {
            double lat = (pt - 0.5*(p0+p1)).dot(perp);
            if (lat >  0.2) sib_above++;
            if (lat < -0.2) sib_below++;
        }
    }

    // Score each side: prefer high SDF, penalise sibling occupation
    double score_above = d_above - 0.3 * sib_above;
    double score_below = d_below - 0.3 * sib_below;

    return (score_above >= score_below) ? +1 : -1;
}

// Build a smooth 2-segment arch: p0(t0) → apex(apex_tan) → p1(t1)
// apex_tan is perpendicular to the p0→p1 direction (or aligned with t0→t1 average)
static BezierCurve buildArch(const Vec2d& p0, const Vec2d& t0,
                               const Vec2d& apex,
                               const Vec2d& p1, const Vec2d& t1,
                               double alpha = 0.4)
{
    // Tangent at apex: average of incoming (p0→apex) and outgoing (apex→p1),
    // clamped to be perpendicular to the net p0→p1 direction for a smooth arch.
    Vec2d leg0 = (apex - p0).norm() > 1e-8 ? (apex - p0).normalized() : t0.normalized();
    Vec2d leg1 = (p1 - apex).norm() > 1e-8 ? (p1 - apex).normalized() : t1.normalized();
    Vec2d apex_tan = (leg0 + leg1);
    if (apex_tan.norm() < 1e-8) apex_tan = leg0;
    apex_tan.normalize();

    BezierSegment s0 = makeCubicG1(p0,    t0.normalized(), apex, apex_tan,  alpha);
    BezierSegment s1 = makeCubicG1(apex,  apex_tan,        p1,   t1.normalized(), alpha);

    BezierCurve c;
    c.segs.push_back(s0);
    c.segs.push_back(s1);
    return c;
}

// Level-1 entry point
static BezierCurve geometricBypass(const Vec2d& p0, const Vec2d& t0,
                                    const Vec2d& p1, const Vec2d& t1,
                                    const SDFField& sdf,
                                    const Polygon2d& fence,
                                    double clearance = 0.3,
                                    const std::vector<std::vector<Vec2d>>& sibling_polys = {})
{
    auto box = probeObstacleAABB(sdf, p0, p1, clearance * 0.8);
    if (!box.valid) {
        // No obstacle detected → direct single-segment curve
        BezierCurve c;
        c.segs.push_back(makeCubicG1(p0, t0.normalized(), p1, t1.normalized(), 0.4));
        return c;
    }

    int side = chooseSide(sdf, box, p0, p1, clearance, sibling_polys);
    Vec2d apex = computeApex(box, p0, p1, side, clearance + 0.5);

    // Verify apex is obstacle-free; if not, try other side
    auto [d_apex, dummy] = sdf.queryWithGrad(apex);
    if (d_apex < clearance) {
        apex = computeApex(box, p0, p1, -side, clearance + 0.5);
    }

    // Verify arch clears obstacles (sample midpoints of each segment)
    BezierCurve arch = buildArch(p0, t0, apex, p1, t1);
    bool arch_clear = true;
    for (auto& seg : arch.segs) {
        for (int i = 1; i < 15; ++i) {
            auto [d, dum] = sdf.queryWithGrad(seg.evaluate((double)i/15));
            if (d < clearance * 0.5) { arch_clear = false; break; }
        }
        if (!arch_clear) break;
    }

    if (!arch_clear) {
        // Try opposite side
        BezierCurve arch2 = buildArch(p0, t0,
            computeApex(box, p0, p1, -side, clearance + 0.5),
            p1, t1);
        bool arch2_clear = true;
        for (auto& seg : arch2.segs) {
            for (int i = 1; i < 15; ++i) {
                auto [d, dum] = sdf.queryWithGrad(seg.evaluate((double)i/15));
                if (d < clearance * 0.5) { arch2_clear = false; break; }
            }
            if (!arch2_clear) break;
        }
        if (!arch2_clear) return {};  // signal Level-2 fallback
        arch = arch2;
    }

    // Validate no intersection with existing sibling curves
    if (!sibling_polys.empty()) {
        auto arch_pts = arch.sampleByArcLength(20);
        for (auto& sp : sibling_polys) {
            for (int ai = 0; ai+1 < (int)arch_pts.size(); ++ai) {
                for (int si = 0; si+1 < (int)sp.size(); ++si) {
                    Vec2d isect;
                    if (segmentsIntersect_internal(
                            arch_pts[ai], arch_pts[ai+1],
                            sp[si],       sp[si+1], &isect)) {
                        // Skip endpoint proximity
                        double de = std::min(
                            (isect-arch_pts.front()).norm(),
                            (isect-arch_pts.back()).norm());
                        if (de > 0.05) { return {}; }  // Level-2 fallback
                            }
                }
            }
        }
    }

    return arch;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Level-2: Full control-point geometry initialisation (方案C 初始化)
//
//  Creates a plausible initial arch WITHOUT grid routing.
//  All control points (including the join point) are free for optimisation.
//  Strategy: tangent-extension intersection gives the arch apex analytically.
// ─────────────────────────────────────────────────────────────────────────────
static BezierCurve geometricInitLevel2(const Vec2d& p0, const Vec2d& t0,
                                        const Vec2d& p1, const Vec2d& t1,
                                        const SDFField& sdf)
{
    // Step 1: find a safe lateral direction (away from obstacle cloud)
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) {
        BezierCurve c;
        c.segs.push_back(makeCubicG1(p0,t0.normalized(),p1,t1.normalized(),0.4));
        return c;
    }
    along = along * (1.0/len);
    Vec2d perp{-along.y(), along.x()};

    // Sample SDF on left vs right of midpath to pick side with more clearance
    double d_left=0, d_right=0;
    int n_probe = 10;
    for (int i=1; i<=n_probe; ++i) {
        double s = (double)i/(n_probe+1)*len;
        Vec2d mid = p0 + s*along;
        auto [dl,dum1] = sdf.queryWithGrad(mid + 1.5*perp);
        auto [dr,dum2] = sdf.queryWithGrad(mid - 1.5*perp);
        d_left += dl; d_right += dr;
    }
    double side = (d_left >= d_right) ? 1.0 : -1.0;

    // Step 2: place 3-waypoint arch
    //   apex at 50% longitudinal, lateral offset = obstacle_radius + clearance
    //   We estimate obstacle_radius from min SDF along direct line
    double min_d = 1e18;
    for (int i=1; i<20; ++i) {
        double t=(double)i/20;
        auto [d,dum] = sdf.queryWithGrad((1-t)*p0+t*p1);
        min_d = std::min(min_d, d);
    }
    double lateral_needed = std::max(1.0, -min_d + 1.5);  // how far to bypass

    Vec2d apex = p0 + 0.5*(p1-p0) + side*lateral_needed*perp;

    // Step 3: build 4-segment arch through: p0 → q1 → apex → q2 → p1
    //   q1 = 25% along, lifted half-way to apex
    //   q2 = 75% along, lifted half-way to apex
    Vec2d q1 = p0 + 0.25*(p1-p0) + side*(lateral_needed*0.6)*perp;
    Vec2d q2 = p0 + 0.75*(p1-p0) + side*(lateral_needed*0.6)*perp;

    // Tangents: Catmull-Rom from waypoints
    std::vector<Vec2d> pts  = {p0, q1, apex, q2, p1};
    std::vector<Vec2d> tans = {t0.normalized(), {}, {}, {}, t1.normalized()};
    for (int i=1; i<=3; ++i) {
        Vec2d d = 0.5*(pts[i+1]-pts[i-1]);
        tans[i] = d.norm()>1e-10 ? d.normalized()
                                  : (pts[i+1]-pts[i]).normalized();
        // Clamp against net direction to avoid loops
        tans[i] = clampTangent(tans[i], along, MAX_TAN_DEV);
    }

    return makeCurveFromKnots(pts, tans, 0.35);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Public API
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve buildInitialCurve(const Vec2d& p0, const Vec2d& t0,
                               const Vec2d& p1, const Vec2d& t1,
                               const SDFField& sdf,
                               const Polygon2d& fence,
                               const std::vector<std::vector<Vec2d>>& sibling_polys)
{
    // U-turn: special case
    if (angleBetween(t0, t1) > M_PI * 0.85)
        return buildTwoSegmentUTurn(p0, t0, p1, t1, sdf, fence);

    // ── Fast path: straight line clear ──────────────────────────────────────
    if (straightLineClear(sdf, p0, p1, 0.15, 30)) {
        BezierCurve c;
        c.segs.push_back(makeCubicG1(p0, t0.normalized(), p1, t1.normalized(), 0.4));
        return c;
    }

    // ── Level-1: geometric direct construction ───────────────────────────────
    {
        BezierCurve arch = geometricBypass(p0, t0, p1, t1, sdf, fence, 0.3, sibling_polys);
        if (!arch.empty()) return arch;
    }

    // ── Level-2: geometry-initialised full-control-point curve ───────────────
    return geometricInitLevel2(p0, t0, p1, t1, sdf);
}

// ─────────────────────────────────────────────────────────────────────────────
//  U-turn (unchanged)
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve buildTwoSegmentUTurn(const Vec2d& p0, const Vec2d& t0,
                                  const Vec2d& p1, const Vec2d& t1,
                                  const SDFField& sdf,
                                  const Polygon2d&)
{
    double d = (p1-p0).norm();
    Vec2d perp{-t0.y(), t0.x()};
    Vec2d mid = 0.5*(p0+p1) + perp*(d*0.5);

    auto [ds, dum1] = sdf.queryWithGrad(mid);
    if (ds < 0.2) {
        Vec2d m2 = 0.5*(p0+p1) - perp*(d*0.5);
        auto [d2, dum2] = sdf.queryWithGrad(m2);
        if (d2 > ds) mid = m2;
    }
    Vec2d tm = (p1-p0).norm()>1e-6 ? (p1-p0).normalized() : t0;
    return makeCurveFromKnots({p0,mid,p1},
                               {t0.normalized(),tm,t1.normalized()}, 0.4);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Legacy: sdfMaxClearancePath retained for any external callers
// ─────────────────────────────────────────────────────────────────────────────
std::vector<Vec2d> sdfMaxClearancePath(const SDFField& sdf,
                                        const Polygon2d& fence,
                                        const Vec2d& start,
                                        const Vec2d& goal,
                                        double /*alpha*/)
{
    // Simplified: just return start→goal; buildInitialCurve no longer uses this.
    (void)sdf; (void)fence;
    return {start, goal};
}
