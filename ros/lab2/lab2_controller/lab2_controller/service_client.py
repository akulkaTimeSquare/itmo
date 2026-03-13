from lab2_interfaces.srv import A368540
from rclpy.node import Node
import rclpy
import time


class ServiceClient(Node):
    def __init__(self, points):
        super().__init__("service_client")
        self.client = self.create_client(A368540, "a368540")
        while not self.client.wait_for_service(timeout_sec=1.0):
            self.get_logger().info('Server not available, waiting...')
        self.points = points
        self.n = len(points)
        self.i = 0
        self.request = A368540.Request()
    
    def send_goal(self):
        if self.i < self.n:
            self.request.x = self.points[self.i][0]
            self.request.y = self.points[self.i][1]
            self.future = self.client.call_async(self.request)


def main(args=None):
    rclpy.init(args=args)
    points = [(2,3), (4,2), (7,1), (8,4), (4,2)]
    service_client = ServiceClient(points=points)

    service_client.get_logger().info(f"Starting our mission...")
    service_client.send_goal()
    while rclpy.ok() and service_client.i < service_client.n:
        rclpy.spin_once(service_client, timeout_sec=0.25)
        if service_client.future.done():
            try:
                response = service_client.future.result()
                if response.success:
                    service_client.get_logger().info(f"Point {points[service_client.i]} reached!")
                    service_client.i += 1
                service_client.send_goal()
            except Exception as e:
                service_client.get_logger().info(f"Service call failed: {e}")
    service_client.get_logger().info(f"All points reached!")
    service_client.destroy_node()
    rclpy.shutdown()


if __name__ == "__main__":
    main()
