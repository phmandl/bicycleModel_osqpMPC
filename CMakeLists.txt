cmake_minimum_required(VERSION 3.10)

# set the project name
project(bicycleModel_osqpMPC)

# add the executable
add_executable(${PROJECT_NAME} bicycleModelMPC.cpp)

# EIGEN STUFF
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

# Find OSQP library and headers
find_package(OsqpEigen REQUIRED)
target_link_libraries(${PROJECT_NAME} OsqpEigen::OsqpEigen)