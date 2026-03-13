import rclpy
from rclpy.node import Node
from sensor_msgs.msg import LaserScan, Image
from geometry_msgs.msg import Twist
from cv_bridge import CvBridge
import cv2
import numpy as np


class DataProcessingNode(Node):
    def __init__(self):
        super().__init__('data_processor')
        
        # parameters
        self.safe_distance = 0.75
        self.final_distance = 0.3
        self.detection_range = 3
        self.turn_duration = 3.0
        self.red_pixel_threshold = 0.15

        # pubs and subs
        self.publisher = self.create_publisher(
            Twist, '/cmd_vel', 10
        )
        self.lidar_sub = self.create_subscription(
            LaserScan, '/diff_drive/scan', self.lidar_callback, 100
        )
        self.camera_sub = self.create_subscription(
            Image, '/camera/image_raw', self.camera_callback, 10
        )

        # for moving
        self.timer = self.create_timer(0.1, self.control_loop)
        self.mode = "forward"
        self.turn_start_time = 0

        # data processing
        self.lidar_data = None
        self.cv_image = None
        self.bridge = CvBridge()

    def lidar_callback(self, msg):
        self.lidar_data = msg

    def camera_callback(self, msg):
        try:
            self.cv_image = self.bridge.imgmsg_to_cv2(msg, "bgr8")
        except Exception as e:
            self.get_logger().error(f"Error converting image: {e}")

    def is_red_detected(self):
        if self.cv_image is None: return 0
        
        hsv = cv2.cvtColor(self.cv_image, cv2.COLOR_BGR2HSV)
        
        mask1 = cv2.inRange(hsv, np.array([0, 120, 70]), np.array([10, 255, 255]))
        mask2 = cv2.inRange(hsv, np.array([170, 120, 70]), np.array([180, 255, 255]))
        
        red_pixels = np.sum(mask1 + mask2 >= 255)
        total_pixels = self.cv_image.shape[0] * self.cv_image.shape[1]
        rt = red_pixels / total_pixels
        
        self.get_logger().info(f'Red ratio: {rt:.2f}')
        return rt > self.red_pixel_threshold

    def get_target_angle(self):
        if self.lidar_data is None: return
        
        ranges = [min(r, 100) for r in self.lidar_data.ranges]
        
        min_dist = min(ranges)
        min_index = ranges.index(min_dist)
        return self.lidar_data.angle_min + min_index * self.lidar_data.angle_increment, min_dist

    def control_loop(self):
        if self.lidar_data is None: return
        
        twist = Twist()
        
        # closest info
        target_angle, min_dist = self.get_target_angle()
        
        if self.mode == "forward":
            twist.linear.x = 1.0
            twist.angular.z = 0.0
            if min_dist <= self.detection_range:
                self.mode = "close"

        elif self.mode == "close":
            if min_dist > self.safe_distance:
                twist.linear.x = 0.5
                twist.angular.z = 4 * target_angle
            else:
                twist.linear.x = 0.0
                twist.angular.z = 0.0
                self.mode = "check"

        elif self.mode == "check":
            if self.is_red_detected():
                self.mode = "final"
            else:
                self.mode = "next"
                self.turn_start_time = self.get_clock().now()

        elif self.mode == "next":
            twist.linear.x = 0.0
            twist.angular.z = min(-0.5, -2.75 * (0.75 - abs(target_angle)))
            
            duration = self.get_clock().now() - self.turn_start_time
            seconds = duration.nanoseconds / 10**9
            
            if seconds > self.turn_duration:
                self.mode = "forward"

        elif self.mode == "final":
            if min_dist > self.final_distance:
                twist.linear.x = 0.5
                twist.angular.z = 3.0 * target_angle 
            else:
                twist.linear.x = 0.0
                twist.angular.z = 0.0
            
        self.publisher.publish(twist)


def main(args=None):
    rclpy.init(args=args)
    node = DataProcessingNode()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
