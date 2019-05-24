#include <ros/ros.h>

#include <forward_image_receiver.h>

int main(int argc, char** argv)
{
    ros::init(argc, argv, "patch_estimator_node");

    ForwardImageReceiver forward_image_receiver;

    // ros::AsyncSpinner spinner(4);
    // spinner.start();
    // ros::waitForShutdown();

    ros::MultiThreadedSpinner spinner(4);
    spinner.spin();

    // ros::spin();
    return 0;
}
