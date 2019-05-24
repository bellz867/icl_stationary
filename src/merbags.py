#!/usr/bin/env python

import rospy
import roslib
import rosbag

bagNew = rosbag.Bag('loop2.bag','w')
bagncr = rosbag.Bag('loop2_ncr.bag')

tripodImage = []
tripodCameraInfo = []
tripodOdomEKF = []
boardOdomEKF = []
turtleOdomEKF = []
turtleOdom = []
turtleImage = []
turtleCameraInfo = []

# get the messages from the ncr bag
for topic,msg,time in bagncr.read_messages(topics=['/camera/image_raw', '/camera/camera_info', '/tripodcam/odomEKF', '/board/odomEKF']):
	#print ('ncr time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
	if topic == '/camera/image_raw':
		# print ('tripod image time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
		tripodImage.append(msg)
	if topic == '/camera/camera_info':
		# print ('tripod camera_info time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
		tripodCameraInfo.append(msg)
	if topic == '/tripodcam/odomEKF':
		# print ('tripod odomEKF time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
		tripodOdomEKF.append(msg)
	if topic == '/board/odomEKF':
		# print ('board odomEKF time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
		boardOdomEKF.append(msg)
	print ('ncr time: ' + '%.8f' % time.to_sec())

bagncr.close()

# get the messages from the turtle bag, wait for the first odomEKF so can adjust the times of the image, camera info. and odom using it
bagturtle = rosbag.Bag('loop2_turtle.bag')
firstOdomEKF = True
adjustTime = rospy.Duration()

for topic,msg,time in bagturtle.read_messages(topics=['/turtlebot3/odom', '/turtlebot3/odomEKF', '/turtlebot3/camera/image_raw', '/turtlebot3/camera/camera_info']):
	#print ('turtle time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
	if (not firstOdomEKF) and topic == '/turtlebot3/odom':
		msg.header.stamp += adjustTime
		# print ('turtle odom time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
		turtleOdom.append(msg)
	if topic == '/turtlebot3/odomEKF':
		# print ('turtle odomEKF time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())

		#if its the first odomEKF subtract the time from stamp time to get adjustment
		if firstOdomEKF:
			adjustTime = msg.header.stamp-time
			firstOdomEKF = False
		turtleOdomEKF.append(msg)

	if (not firstOdomEKF) and topic == '/turtlebot3/camera/image_raw':
		msg.header.stamp += adjustTime
		# print ('turtle image time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
		turtleImage.append(msg)
	if (not firstOdomEKF) and topic == '/turtlebot3/camera/camera_info':
		msg.header.stamp += adjustTime
		# print ('turtle camera_info time: ' + '%.8f' % time.to_sec() + ' topic: ' + topic + ' type: ' + msg._type + ' stamp: ' + '%.8f' % msg.header.stamp.to_sec())
		turtleCameraInfo.append(msg)
	print ('turtle time: ' + '%.8f' % time.to_sec())
bagturtle.close()

bagTopics = ['/camera/image_raw', '/camera/camera_info', '/tripodcam/odomEKF', '/board/odomEKF', '/turtlebot3/odom', '/turtlebot3/odomEKF', '/turtlebot3/camera/image_raw', '/turtlebot3/camera/camera_info']
bagMsgs = [tripodImage,tripodCameraInfo,tripodOdomEKF,boardOdomEKF,turtleOdom,turtleOdomEKF,turtleImage,turtleCameraInfo]

# go through the lists and while there is still something in always grab the oldest message and add it to the new bag
while len(bagMsgs) > 0:
	minTime = rospy.Time()
	firstMsg = True
	minTimeInd = 0
	TimeInd = 0
	removeInds = []

	while TimeInd < len(bagTopics):
		if len(bagMsgs[TimeInd]) > 0:
			msg = bagMsgs[TimeInd][0]
			if firstMsg:
				minTime = msg.header.stamp
				firstMsg = False
				minTimeInd = TimeInd
			else:
				Time = msg.header.stamp
				if (minTime-Time).to_sec() > 0.0:
					minTimeInd = TimeInd
					minTime = Time
		else:
			removeInds.append(TimeInd)
		TimeInd += 1

	try:
		print ('save time: ' + '%.8f' % bagMsgs[minTimeInd][0].header.stamp.to_sec() + ' topic: ' + bagTopics[minTimeInd])

		if len(bagMsgs[minTimeInd]) > 0:
			# take the minimum and put it in the bag
			bagNew.write(bagTopics[minTimeInd],bagMsgs[minTimeInd].pop(0),minTime)
	except:
		bagNew.close()


	for removeInd in removeInds:
		bagTopics.pop(removeInd)
		bagMsgs.pop(removeInd)

bagNew.close()
