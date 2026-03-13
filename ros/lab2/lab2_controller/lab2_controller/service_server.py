from lab2_interfaces.srv import A368540
from rclpy.node import Node
import rclpy
from geometry_msgs.msg import Twist
from turtlesim.msg import Pose
import math


class ServiceServer(Node):
    def __init__(self, distance_tolerance=0.01, angle_tolerance=0.05):
        super().__init__("service_server")
        self.subscription = self.create_subscription(
            Pose,
            "turtle1/pose",
            self.pose_callback,
            10
        )
        self.current_pose = None
        self.publisher = self.create_publisher(Twist, "turtle1/cmd_vel", 10)
        self.server = self.create_service(A368540, "a368540", self.execute_callback)
        self.distance_tolerance = distance_tolerance
        self.angle_tolerance = angle_tolerance
    
    def pose_callback(self, pose: Pose):
        self.current_pose = pose

    def execute_callback(self, request, response):
        if self.current_pose is None:
            response.success = False
            rclpy.spin_once(self, timeout_sec=0.05)
            return response
            
        target_x = request.x
        target_y = request.y

        distance_error = ((self.current_pose.x - target_x) ** 2 + (self.current_pose.y - target_y) ** 2) ** 0.5
        if distance_error < self.distance_tolerance:
            response.success = True
            self.get_logger().info(f"Point ({target_x}, {target_y}) reached!")
            return response
        
        cmd = self.createTwist(target_x, target_y)
        self.publisher.publish(cmd)
        
        rclpy.spin_once(self, timeout_sec=0.15)
        distance_updated = ((self.current_pose.x - target_x) ** 2 + (self.current_pose.y - target_y) ** 2) ** 0.5
        if distance_updated < self.distance_tolerance:
            response.success = True
            self.get_logger().info(f"Point ({target_x}, {target_y}) reached!")
        else:
            response.success = False
        return response

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
    service_server = ServiceServer()
    rclpy.spin(service_server)
    rclpy.shutdown()


if __name__ == "__main__":
    main()
