cmake_minimum_required(VERSION 3.8)
project(uclv_systems_lib)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)

find_package (Eigen3 REQUIRED NO_MODULE)

include_directories(include)

# Export old-style CMake variables
ament_export_include_directories("include")
# ament_export_libraries(${PROJECT_NAME})

# Export modern CMake targets
# ament_export_targets(${PROJECT_NAME}Targets HAS_LIBRARY_TARGET)

# specific order: dependents before dependencies
ament_export_dependencies(Eigen3)

install(
  DIRECTORY include/
  DESTINATION include
)

# install(
#   TARGETS ${PROJECT_NAME} 
#   EXPORT ${PROJECT_NAME}Targets
#   LIBRARY DESTINATION lib
#   ARCHIVE DESTINATION lib
#   RUNTIME DESTINATION bin
#   INCLUDES DESTINATION include
# )

# ament_generate_version_header(${PROJECT_NAME})

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
  # add_executable(test_compile test/test_compile.cpp)
  # target_link_libraries(test_compile Eigen3::Eigen)
endif()

ament_package()
