cmake_minimum_required(VERSION 2.8.3)
project(icl_stationary)

set(CMAKE_CXX_FLAGS "-std=c++11 -msse -msse2 -msse3 ${CMAKE_CXX_FLAGS}")

find_package(cmake_modules REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  tf
  cv_bridge
  image_transport
  image_geometry
  sensor_msgs
  geometry_msgs
  message_generation
  std_msgs
  nav_msgs
  pcl_conversions
  pcl_ros
)

add_message_files(
  FILES
  Output.msg
  Roi.msg
)

generate_messages(
  DEPENDENCIES
  std_msgs
  geometry_msgs
  nav_msgs
)

catkin_package(
  CATKIN_DEPENDS message_runtime
)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  "include"
)

add_executable(patch_estimator_node_station src/depth_estimator.cpp src/depth_estimator_ekf.cpp src/depth_estimator_ls.cpp src/depth_estimator_norm.cpp src/depth_estimator_icl.cpp src/depth_estimator_icl_ext.cpp src/vector_derivative_estimator.cpp src/feature_derivative_estimator.cpp src/patch_estimator.cpp src/data_save.cpp src/forward_image_receiver.cpp src/helper_functions.cpp src/patch_estimator_main.cpp)
#add_executable(wall_mapper_node src/keyframe_planes.cpp src/wall_mapper.cpp src/wall_mapper_main.cpp)
add_executable(odom_estimator_node_station src/odom_estimator.cpp src/odom_estimator_main.cpp src/helper_functions.cpp)

add_dependencies(patch_estimator_node_station ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
#add_dependencies(wall_mapper_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
add_dependencies(odom_estimator_node_station ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

target_link_libraries(patch_estimator_node_station
  ${catkin_LIBRARIES}
  ${OpenCV_LIBRARIES}
)

#target_link_libraries(wall_mapper_node
#  ${catkin_LIBRARIES}
#  ${OpenCV_LIBRARIES}
#)

target_link_libraries(odom_estimator_node_station
  ${catkin_LIBRARIES}
  ${OpenCV_LIBRARIES}
)
