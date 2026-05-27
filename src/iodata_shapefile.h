//
// Created by ubuntu on 26-5-24.
//

#ifndef IODATA_SHAPEFILE_H
#define IODATA_SHAPEFILE_H
#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <filesystem>
#include "shapefile.hpp"
#include "types.h"

inline std::vector<ShapePoint> toArray(const std::vector<Vec2d>& points) {
    std::vector<ShapePoint> shppoints;
    for (auto p : points) {
        shppoints.emplace_back(ShapePoint{p.x(), p.y(), 0.0});
    }
    return shppoints;
};

static bool save(IntersectionInput& inp, std::string dir, std::string prefix) {
    auto writeLanes = [&](const std::string& dir, const std::string& fname,
        const std::vector<Lane>& lanes) {
        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"LANE_ORDER",  'C', 64, 0},
            {"GROUP_ID",'C', 64, 0},
            {"GROUP_TYPE",'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < lanes.size(); ++i) {
            const auto& pl = lanes[i];
            std::vector<std::string> attrs = {pl.id, "-1", "", ""};
            std::vector<ShapePoint> points = toArray(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeLaneEdges = [&](const std::string& dir, const std::string& fname,
        const std::vector<LaneEdge>& edgelines){
        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"GROUP_ID",'C', 64, 0},
            {"GROUP_TYPE",'C', 64, 0},
            {"LEFT_CLINE_ID",'C', 64, 0},
            {"RIGHT_CLINE_ID",'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < edgelines.size(); ++i) {
            const auto& pl = edgelines[i];
            std::vector<std::string> attrs = {pl.id, "-1", "", "", ""};
            std::vector<ShapePoint> points = toArray(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeStopLines = [&](const std::string& dir, const std::string& fname,
        const std::vector<StopLine>& stoplines){
        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"ENTRY_GROUP_ID",'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < stoplines.size(); ++i) {
            const auto& pl = stoplines[i];
            std::vector<std::string> attrs = {pl.id, pl.associated_group_id};
            std::vector<ShapePoint> points = toArray(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeRoadEdges = [&](const std::string& dir, const std::string& fname,
        const std::vector<Boundary>& boundaries){
        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"TYPE",'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < boundaries.size(); ++i) {
            const auto& pl = boundaries[i];
            std::vector<std::string> attrs = {pl.id, std::to_string((int)pl.type)};
            std::vector<ShapePoint> points = toArray(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeObstacles = [&](const std::string& dir, const std::string& fname,
        const std::vector<Obstacle>& obstacles){
        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"TYPE",'C', 64, 0},
        };
        int shpType = SHP_POLYGONZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < obstacles.size(); ++i) {
            const auto& pg = obstacles[i];
            std::vector<std::string> attrs = {pg.id, "-1"};//std::to_string(pg.type)
            std::vector<ShapePoint> points = toArray(pg.geometry.outer);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeCrosswalks = [&](const std::string& dir, const std::string& fname,
        const std::vector<Crosswalk>& crosswalks){
        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"TYPE",'C', 64, 0},
        };
        int shpType = SHP_POLYGONZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < crosswalks.size(); ++i) {
            const auto& pg = crosswalks[i];
            std::vector<std::string> attrs = {pg.id, "-1"};//std::to_string(pg.type)
            std::vector<ShapePoint> points = toArray(pg.geometry.outer);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    writeLanes(dir, prefix + "_lanes", inp.lanes);
    writeLaneEdges(dir, prefix + "_laneedges", inp.lane_edges);
    writeStopLines(dir, prefix + "_stoplines", inp.stop_lines);
    writeRoadEdges(dir, prefix + "_roadedges", inp.boundaries);
    writeObstacles(dir, prefix + "_obstacles", inp.obstacles);
    writeCrosswalks(dir, prefix + "_crosswalks", inp.crosswalks);
    return true;
}

static bool save(IntersectionOutput& out, std::string dir, std::string prefix) {
    auto writeLanes = [&](const std::string& dir, const std::string& fname,
        const std::vector<ConnectivityCurve>& lanes) {
        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"TURN_TYPE",  'C', 64, 0},
            {"FLANE",'C', 64, 0},
            {"TLANE",'C', 64, 0},
        };
        int gid = -1;
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < lanes.size(); ++i) {
            const auto& pl = lanes[i];
            if (pl.curve->numSegments() == 0) continue;
            std::vector<std::string> attrs = {pl.id, std::to_string((int)pl.turn_type), pl.entry_lane_id, pl.exit_lane_id};
            std::vector<ShapePoint> points = toArray(pl.curve->sample(50));
            ShapeRecord record = {
                ++gid, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeLaneEdges = [&](const std::string& dir, const std::string& fname,
        const std::vector<LaneEdge>& edgelines){
        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"GROUP_ID",'C', 64, 0},
            {"GROUP_TYPE",'C', 64, 0},
            {"LEFT_CLINE_ID",'C', 64, 0},
            {"RIGHT_CLINE_ID",'C', 64, 0},
        };
        int shpType = SHP_POLYLINEZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < edgelines.size(); ++i) {
            const auto& pl = edgelines[i];
            std::vector<std::string> attrs = {pl.id, "-1", "", "", ""};
            std::vector<ShapePoint> points = toArray(pl.geometry.points);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    auto writeAreas = [&](const std::string& dir, const std::string& fname,
        const IntersectionArea& inter_area){
        std::vector<Polygon2d> areas = {inter_area.fine_area, inter_area.coarse_area};

        std::filesystem::create_directories(dir);
        std::vector<DbfField> fields = {
            {"ID",       'C', 64, 0},
            {"TYPE",'C', 64, 0},
        };
        int shpType = SHP_POLYGONZ;
        std::vector<ShapeRecord> records = {};
        for (int i = 0; i < areas.size(); ++i) {
            const auto& pg = areas[i];
            std::vector<std::string> attrs = {std::to_string(i), "-1"};//std::to_string(pg.type)
            std::vector<ShapePoint> points = toArray(pg.outer);
            ShapeRecord record = {
                i, shpType, attrs, points, {0}
            };
            records.emplace_back(record);
        }
        ShapefileEngine::write(dir, fname, shpType, fields, records);
        //std::cout << "[INFO]  Shapefile written: " << dir << "/" << fname << ".*\n";
    };
    writeLanes(dir, prefix + "_lanes", out.connectivity_curves);
    writeLaneEdges  (dir, prefix + "_laneedges", out.lane_edges);
    writeAreas    (dir, prefix + "_areas", out.area);
#ifdef PROJECT_ROOT_DIR
    std::string proj_dir = PROJECT_ROOT_DIR;
    std::cout << "PROJECT_ROOT_DIR : " << proj_dir << std::endl;
    std::string qgis_fname = "/intersection_gen.qgs";
    std::string qgis_temp = proj_dir + qgis_fname;
    std::string dest_file = std::filesystem::path(dir).parent_path().parent_path().string() + qgis_fname;
    if (std::filesystem::exists(qgis_temp))
        std::filesystem::copy_file(qgis_temp, dest_file,
            std::filesystem::copy_options::skip_existing);
#endif
    return true;
}
#endif //IODATA_SHAPEFILE_H
