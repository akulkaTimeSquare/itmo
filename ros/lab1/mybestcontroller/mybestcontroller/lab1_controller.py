#!/usr/bin/env python3
import rclpy
from rclpy.node import Node
from turtlesim.msg import Pose
from geometry_msgs.msg import Twist
from builtin_interfaces.msg import Duration
from math import atan2, pi


class ControllerNode(Node):
    def __init__(self, points, distance_tolerance=0.01, angle_tolerance=0.05):
        super().__init__('controller_node')
        if self.get_namespace().startswith("/ns1"):
            self.subscription = self.create_subscription(
                Pose,
                'turtle1/pose',
                self.pose_callback_turtle1,
                10)
            self.publisher_time = self.create_publisher(Duration, 'result_368540', 10)
            
            self.distance_tolerance = distance_tolerance
            self.angle_tolerance = angle_tolerance
            
            self.created_time = self.get_clock().now()
            self.points = points
            self.n = len(points)
            self.i = 0

        if self.get_namespace().startswith("/ns2"):
            ns = self.get_namespace()
            ns = ns[:3] + "1" + ns[4:]
            self.get_logger().info(ns)
            self.subscription = self.create_subscription(
                Pose,
                ns + '/turtle1/pose',
                self.pose_callback_turtle2,
                10)

        self.publisher = self.create_publisher(Twist, 'turtle1/cmd_vel', 10)
        

    def make_cmd_to_move(self, pose: Pose, destination_point):
        cmd = Twist()

        dx = destination_point[0] - pose.x
        dy = destination_point[1] - pose.y
        angle_to_goal = atan2(dy, dx)

        distance = (dx ** 2 + dy ** 2) ** 0.5
        angle_diff = angle_to_goal - pose.theta
        angle_diff = (angle_diff + pi) % (2 * pi) - pi

        if distance < self.distance_tolerance:
            return None
          
        if abs(angle_diff) > self.angle_tolerance:
            cmd.angular.z = min(angle_diff, 1.0) if angle_diff > 0 else max(angle_diff, -1.0)
            cmd.linear.x = 0.0
        else:
            cmd.linear.x = min(distance, 1.0)
            cmd.angular.z = angle_diff

        return cmd

    def pose_callback_turtle1(self, pose: Pose):
        if self.i >= self.n:
            return
        cmd = self.make_cmd_to_move(pose, self.points[self.i])
        if cmd is None:
            self.i += 1
            time_diff = self.get_clock().now() - self.created_time
            time_msg = time_diff.to_msg()
            if self.i >= self.n:
                self.publisher_time.publish(time_msg)
                self.get_logger().info(f"All points reached!")
            else:
                self.get_logger().info(f"Reached point {self.points[self.i-1]}, moving to next point...")
            self.get_logger().info(f"Time: {time_msg.sec} sec and {time_msg.nanosec} nanosec")
        else:
            self.publisher.publish(cmd)
    
    def pose_callback_turtle2(self, pose: Pose):
        cmd = Twist()
        cmd.linear.x = pose.linear_velocity 
        cmd.angular.z = pose.angular_velocity
        self.publisher.publish(cmd)


def main(args=None):
    rclpy.init(args=args)
    node = ControllerNode(points=[(2,3), (4,2), (7,1), (8,4), (4,2)])
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        if rclpy.ok():
            rclpy.shutdown()


if __name__ == '__main__':
    main()