add_library(intersection_gen_lib STATIC
    intersection_shape_generator.cpp
        clipper.cpp

    curve/bezier.cpp
    curve/hermite_init.cpp
    curve/curve_utils.cpp

    optimizer/sdf_field.cpp
    optimizer/penalty_cost.cpp
    optimizer/lbfgs_solver.cpp

    constraints/infeasibility_detector.cpp
    constraints/cluster_order.cpp

        generator/connectivity_generator.cpp
        generator/lane_edge_generator.cpp
        generator/fine_area_generator.cpp
        generator/polygon_builder.h
        generator/edge_line_generator.h
)

target_include_directories(intersection_gen_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(DEPENDS_EXISTS STREQUAL "TRUE")
    target_link_libraries(intersection_gen_lib PUBLIC)
else ()
    target_link_libraries(intersection_gen_lib PUBLIC Eigen3::Eigen)
endif ()

target_compile_options(intersection_gen_lib PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -O2>
    $<$<CXX_COMPILER_ID:MSVC>:/W3 /O2>
)
