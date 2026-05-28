#pragma once
#include <Eigen/Dense>
#include <optional>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <functional>
#include <cmath>
#include <algorithm>
#include <limits>

using Vec2d = Eigen::Vector2d;
using VecXd = Eigen::VectorXd;

// ── Geometry primitives ──────────────────────────────────────
struct BoundingBox2d {
    Vec2d min_pt{ 1e18, 1e18};
    Vec2d max_pt{-1e18,-1e18};
    void  expand(const Vec2d& p){min_pt=min_pt.cwiseMin(p);max_pt=max_pt.cwiseMax(p);}
    bool  intersects(const BoundingBox2d& o)const{
        return max_pt[0]>=o.min_pt[0]&&min_pt[0]<=o.max_pt[0]&&
               max_pt[1]>=o.min_pt[1]&&min_pt[1]<=o.max_pt[1];}
    bool  contains(const Vec2d& p)const{
        return p[0]>=min_pt[0]&&p[0]<=max_pt[0]&&
               p[1]>=min_pt[1]&&p[1]<=max_pt[1];}
    double width() const{return max_pt[0]-min_pt[0];}
    double height()const{return max_pt[1]-min_pt[1];}
    bool   empty() const{return min_pt[0]>max_pt[0];}
};

struct LineString2d {
    std::vector<Vec2d> points;
    BoundingBox2d bbox()const{BoundingBox2d b;for(auto&p:points)b.expand(p);return b;}
};

struct Polygon2d {
    std::vector<Vec2d>              outer;
    std::vector<std::vector<Vec2d>> holes;
    bool empty()const{return outer.empty();}
    BoundingBox2d bbox()const{BoundingBox2d b;for(auto&p:outer)b.expand(p);return b;}
    bool contains(const Vec2d& p)const; // defined in fence_check.cpp
};

// ── Bézier types (defined here to avoid circular deps) ───────
struct BezierSegment {
    std::array<Vec2d, 4> ctrl;
    // B(t) 求值
    Vec2d  evaluate(double t)const;

    Vec2d  tangent (double t)const;

    // 曲率 κ(t) = |B'×B''| / |B'|³
    double curvature(double t)const;

    std::pair<BezierSegment,BezierSegment> splitAt(double t)const;

    BoundingBox2d bbox()const;

    // 弧长（数值积分 - 辛普森法）
    double arcLength(int samples=20)const;

    double arcLengthToParam(double s,int samples=50)const;

    // B'(t) 一阶导数
    Vec2d evalDeriv1(double t) const;

    // B''(t) 二阶导数
    Vec2d evalDeriv2(double t) const;

    // 最大曲率（采样估算）
    double maxCurvature(int samples = 50) const;

    // 固定数量采样
    std::vector<Vec2d> sampleCount(int n) const;

    // 固定间距采样（弧长近似）
    std::vector<Vec2d> sampleBySpacing(double spacing) const;

    // 自适应采样（基于角度误差）
    std::vector<Vec2d> sampleAdaptive(double maxAngleDeg, double maxSeg, double minSeg) const;
    void subdivide(double t0, double t1, double maxAngle, double maxSeg, double minSeg,
               std::vector<Vec2d>& pts, int depth) const;

    // 获取 α（P1相对P0的标量，T0方向）
    double getAlpha(const Vec2d& T0) const;

    double getBeta(const Vec2d& T3) const;
};

struct BezierCurve {
    std::vector<BezierSegment> segs;
    bool   empty()const{return segs.empty();}
    int    numSegments()const{return (int)segs.size();}
    Vec2d  startPt() const;
    Vec2d  endPt()   const;
    Vec2d  startTan()const;
    Vec2d  endTan()  const;
    Vec2d  evaluate(double u)const;
    Vec2d  tangent (double u)const;
    std::vector<Vec2d> sample(int n)const;
    std::vector<Vec2d> sampleByArcLength(int n)const;
    double arcLength()const;
    BoundingBox2d bbox()const;
    double maxCurvature(int sps=20)const;
};

// ── Lane IDs ─────────────────────────────────────────────────
using LaneId      = std::string;
using LaneGroupId = std::string;
using LaneEdgeId  = std::string;
using ConnId      = std::string;
using InterId     = std::string;
using AttrMap = std::map<std::string, std::string>;

struct LaneEdge {
    LaneEdgeId   id;
    LineString2d geometry;
    bool         is_shared=false;
    std::optional<std::pair<LaneId,LaneId>> shared_by;

    AttrMap attrs;
    // Vec2d connectionPt;
    // Vec2d tangentDir;
    std::string groupId;
    int lineOrder = 0; // 组内横向排序（0=最内侧）
};

struct Lane {
    LaneId id;
    LaneEdgeId left_edge_id;
    LaneEdgeId right_edge_id;
    double width = 3.5;
    LineString2d geometry;

    AttrMap attrs;
    // Vec2d connectionPt; // 连接点坐标
    // Vec2d tangentDir; // 连接点切线（指向路口内侧）
    std::string groupId;
    int laneOrder = 0; // 组内横向排序（0=最内侧）
};

enum class GroupRole{Entry,Exit};

struct LaneGroup {
    LaneGroupId             id;
    GroupRole               role;
    std::vector<LaneId>     lanes;
    std::vector<LaneEdgeId> boundaries;
    // NOTE: direction & ref_point are retained for legacy reasons but are
    // NOT used by curve generation — all endpoints & tangents come from
    // Lane::geometry directly.
    //Vec2d direction{1,0},ref_point{0,0};

    // added fields;
    AttrMap attrs;
};

enum class TurnType{UTurnLeft=0,TurnLeft,Straight,TurnRight,UTurnRight};

inline const char* turnTypeName(TurnType t){
    switch(t){
    case TurnType::UTurnLeft: return"UTurnLeft";
    case TurnType::TurnLeft:  return"TurnLeft";
    case TurnType::Straight:  return"Straight";
    case TurnType::TurnRight: return"TurnRight";
    case TurnType::UTurnRight:return"UTurnRight";}return"Unknown";}

struct DivergeMergeInfo{
    enum class Type{None,Diverge,Merge}type=Type::None;
    double offset_arc_length=0,lateral_offset=0;
};

struct Connectivity{
    ConnId id;
    LaneId entry_lane_id,exit_lane_id;
    TurnType turn_type=TurnType::Straight;
    DivergeMergeInfo dm_info;

    LaneGroupId enterGroupId;
    LaneGroupId exitGroupId;
};

struct Obstacle {
    std::string id;
    Polygon2d geometry, buffered_geometry;
};

enum class BoundaryType{RoadEdge,MedianStrip,GreenBelt,Other};

struct Boundary{
    std::string id;
    BoundaryType type=BoundaryType::RoadEdge;
    LineString2d geometry;
};

struct StopLine{
    std::string id;
    LineString2d geometry;
    LaneGroupId associated_group_id;
    Vec2d normal_direction{0,1};
};

struct Crosswalk{
    std::string id;
    Polygon2d geometry;
    Vec2d crossing_direction{0,1};
};

// 路口面
struct IntersectionArea {
    InterId id;
    Polygon2d geometry;
    bool is_rough; //true:粗糙路口面，false:精细路口面
};

enum class CurveStatus{OK,WarnA2,Degraded,Infeasible};

struct ViolationInfo{
    enum class InfeasibilityType{None,NarrowPassage,Sandwich,TopologicalBlock,ForcedCross};
    InfeasibilityType type=InfeasibilityType::None;
    double max_obstacle_penetration=0,max_fence_overflow=0,fence_expansion_applied=0;
    std::vector<Vec2d> exempt_crosses;
    std::string reason;
};

struct ConnectivityCurve{
    ConnId id;
    LaneId entry_lane_id;
    LaneId exit_lane_id;
    TurnType turn_type=TurnType::Straight;
    std::optional<BezierCurve> curve;   // BezierCurve is complete above
    CurveStatus status=CurveStatus::OK;
    ViolationInfo violation;

    LaneEdgeId left_edge_id = "";
    LaneEdgeId right_edge_id = "";
};

struct ConnectivityLaneEdge {
    LaneEdgeId id;
    LineString2d geometry;
    bool is_shared = false;
    std::optional<std::pair<LaneId, LaneId>> shared_by;//边线左侧车道，边线右侧车道

    AttrMap attrs;
    // Vec2d connectionPt;
    // Vec2d tangentDir;
    std::string groupId;
    int lineOrder = 0; // 组内横向排序（0=最内侧）
};


struct IntersectionOutput{
    std::vector<ConnectivityCurve>      connectivity_curves;
    std::vector<ConnectivityLaneEdge>   lane_edges;
    IntersectionArea                    area;
    struct PerfStats{
        double sdf_build_ms=0,precheck_ms=0,optimize_ms=0,
               smooth_ms=0,edge_gen_ms=0,area_gen_ms=0;
    }perf;
};

// ── Lane geometry helpers ───────────────────────────────────
// Entry lane: geometry runs outside→intersection; last point is at intersection edge.
// Exit  lane: geometry runs intersection→outside; first point is at intersection edge.

inline Vec2d laneEntryPoint(const LineString2d& geom) {
    if (geom.points.size() < 1) return Vec2d(0,0);
    return geom.points.back();
}
inline Vec2d laneEntryTangent(const LineString2d& geom) {
    int n = (int)geom.points.size();
    if (n < 2) return Vec2d(1,0);
    Vec2d d = geom.points[n-1] - geom.points[n-2];
    return d.norm() > 1e-10 ? d.normalized() : Vec2d(1,0);
}
inline Vec2d laneExitPoint(const LineString2d& geom) {
    if (geom.points.size() < 1) return Vec2d(0,0);
    return geom.points.front();
}
inline Vec2d laneExitTangent(const LineString2d& geom) {
    int n = (int)geom.points.size();
    if (n < 2) return Vec2d(1,0);
    Vec2d d = geom.points[1] - geom.points[0];
    return d.norm() > 1e-10 ? d.normalized() : Vec2d(1,0);
}

struct IntersectionInput{
    std::vector<LaneGroup>    lane_groups;
    std::vector<Lane>         lanes;
    std::vector<LaneEdge>     lane_edges;
    std::vector<Connectivity> connectivities;
    std::vector<Obstacle>     obstacles;
    std::vector<Boundary>     boundaries;
    std::vector<StopLine>     stop_lines;
    std::vector<Crosswalk>    crosswalks;
    IntersectionArea          area;

    std::string id; //路口ID

    const bool IsEntryLaneEdge(const LaneEdgeId& id) const;
    const bool IsEntryLane(const LaneId& id) const;
    const Lane*      findLane (const LaneId&)const;
    const LaneGroup* findGroup(const LaneGroupId&)const;
    const LaneEdge*  findEdge (const LaneEdgeId&)const;
    bool laneGroupExists(const LaneGroupId&)const;
    std::pair<Vec2d,Vec2d> entryPtDir(const LaneId&)const;
    std::pair<Vec2d,Vec2d> exitPtDir (const LaneId&)const;
};

// ── AdaptiveRefineResult forward-declared here ───────────────
struct SDFField; // forward
struct AdaptiveRefineResult {
    BezierCurve curve;
    bool was_split = false;
};
