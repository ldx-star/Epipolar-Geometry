cmake_minimum_required(VERSION 3.5)
project(Epipolar_Geometry)
set(CMAKE_CXX_STANDARD 11)

#add_library(lib test_camera.cpp)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

find_package(OpenCV REQUIRED)
include_directories(${OPENCV_INCLUDE_DIRS})


add_executable(Test_Camera test_camera.cpp)
add_executable(Test_Fundamental_8_point test_fundamental_8_point.cpp)
add_executable(Test_Fundamental_Ransac test_fundamental_ransac.cpp)
add_executable(Test_Pose_Estimation test_pose_estimation.cpp)

target_link_libraries(Test_Camera ${OpenCV_LIBS})
target_link_libraries(Test_Fundamental_8_point ${OpenCV_LIBS})
target_link_libraries(Test_Fundamental_Ransac ${OpenCV_LIBS})
target_link_libraries(Test_Pose_Estimation ${OpenCV_LIBS})
