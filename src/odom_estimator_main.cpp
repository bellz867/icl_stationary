#include <ros/ros.h>

#include <odom_estimator.h>

int main(int argc, char** argv)
{
    ros::init(argc, argv, "odom_estimator_node");

    OdomEstimator odomEstimator;

    //ros::AsyncSpinner spinner(4);
    //spinner.start();
    //ros::waitForShutdown();

    //ros::MultiThreadedSpinner spinner(4);
    //spinner.spin();

    ros::spin();
    return 0;
}
