from lab2_interfaces.action import A368540
from geometry_msgs.msg import Twist
from turtlesim.msg import Pose
from rclpy.action import ActionServer
from rclpy.node import Node
import rclpy
import math


class ActionServerNode(Node):
    def __init__(self, distance_tolerance=0.01, angle_tolerance=0.05):
        super().__init__("action_server")
        self.publisher = self.create_publisher(Twist, "turtle1/cmd_vel", 10)
        self.subscription = self.create_subscription(
            Pose,
            "turtle1/pose",
            self.pose_callback,
            10
        )
        self.server = ActionServer(
            self,
            A368540,
            "a368540",
            self.execute_callback
        )
        self.current_pose = None
        self.distance_tolerance = distance_tolerance
        self.angle_tolerance = angle_tolerance
        self.last_time_published = None

    def pose_callback(self, pose: Pose):
        self.current_pose = pose
    
    def execute_callback(self, goal_handle):
        if self.current_pose is None:
            goal_handle.abort()
            rclpy.spin_once(self, timeout_sec=0.05)
            return A368540.Result(success=False)
        
        target_x = goal_handle.request.x
        target_y = goal_handle.request.y
        feedback = A368540.Feedback()
        distance_error = ((self.current_pose.x - target_x) ** 2 + (self.current_pose.y - target_y) ** 2) ** 0.5
        while distance_error > self.distance_tolerance:
            if goal_handle.is_cancel_requested:
                goal_handle.canceled()
                self.get_logger().info("Goal canceled!")
                return A368540.Result(success=False)
                
            if self.last_time_published is None or \
             self.get_clock().now() - self.last_time_published > rclpy.duration.Duration(seconds=1.5):
                self.last_time_published = self.get_clock().now()
                feedback.distance_error = distance_error
                goal_handle.publish_feedback(feedback)
                self.get_logger().info(f"Progress: distance = {distance_error:.4f}")

            cmd = self.createTwist(target_x, target_y)
            self.publisher.publish(cmd)
            rclpy.spin_once(self, timeout_sec=0.05)
            distance_error = ((self.current_pose.x - target_x) ** 2 + (self.current_pose.y - target_y) ** 2) ** 0.5
        
        goal_handle.succeed()
        self.get_logger().info(f"Point ({target_x}, {target_y}) reached!")
        return A368540.Result(success=True)

    def createTwist(self, target_x, target_y):
        cmd = Twist()
        dx = target_x - self.current_pose.x
        dy = target_y - self.current_pose.y
        angle_to_goal = math.atan2(dy, dx)
        angle_diff = angle_to_goal - self.current_pose.theta
        angle_diff = (angle_diff + math.pi) % (2 * math.pi) - math.pi
        distance = (dx ** 2 + dy ** 2) ** 0.5

        if abs(angle_diff) > self.angle_tolerance:
            cmd.angular.z = min(angle_diff, 1.0) if angle_diff > 0 else max(angle_diff, -1.0)
            cmd.linear.x = 0.0
        else:
            cmd.linear.x = min(distance, 1.0)
            cmd.angular.z = angle_diff
        return cmd


def main(args=None):
    rclpy.init(args=args)
    action_server = ActionServerNode()
    rclpy.spin(action_server)
    rclpy.shutdown()


if __name__ == "__main__":
    main()