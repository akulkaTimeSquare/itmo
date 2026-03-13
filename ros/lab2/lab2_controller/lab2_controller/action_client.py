from lab2_interfaces.action import A368540
from rclpy.action import ActionClient
from rclpy.node import Node
import rclpy


class ActionClientNode(Node):
    def __init__(self, points):
        super().__init__("action_client")
        self.client = ActionClient(self, A368540, "a368540")
        self.points = points
        self.n = len(points)
        self.i = 0
    
    def send_goal(self):
        while not self.client.wait_for_server(timeout_sec=1.0):
            self.get_logger().info('Server not available, waiting...')

        if self.i >= self.n:
            self.get_logger().info(f"All points reached!")
            self.destroy_node()
            rclpy.shutdown()
            return

        target_x = self.points[self.i][0]
        target_y = self.points[self.i][1]

        goal_msg = A368540.Goal()
        goal_msg.x = target_x
        goal_msg.y = target_y
        self.get_logger().info(f"New goal: point({target_x}, {target_y})")
        self.client.send_goal_async(
            goal_msg,
            feedback_callback=self.feedback_callback
        ).add_done_callback(self.goal_response_callback)
    
    def goal_response_callback(self, future):
        goal_handle = future.result()
        if not goal_handle.accepted:
            self.get_logger().warn(f"Goal canceled!")
            return
        
        self._response_handle = goal_handle.get_result_async()
        self._response_handle.add_done_callback(self.result_response_callback)
    
    def result_response_callback(self, future):
        result = future.result().result
        if result.success:
            self.get_logger().info(f"Mission completed!")
            self.i += 1
        else:
            self.get_logger().warn(f"Failed! :(")
        self.send_goal()
    
    def feedback_callback(self, feedback_msg):
        distance_error = feedback_msg.feedback.distance_error
        self.get_logger().info(f"Feedback: distance = {distance_error:.4f}")


def main(args=None):
    rclpy.init(args=args)
    points = [(2,3), (4,2), (7,1), (8,4), (4,2)]
    action_client = ActionClientNode(points=points)
    action_client.send_goal()
    rclpy.spin(action_client)


if __name__ == "__main__":
    main()