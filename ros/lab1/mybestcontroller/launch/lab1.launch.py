from launch import LaunchDescription
from launch_ros.actions import Node

def generate_launch_description():
    return LaunchDescription([
        Node(
            package='turtlesim',
            executable='turtlesim_node',
            name='sim1',
            output='screen',
            namespace='ns1_368540'
        ),
        Node(
            package='mybestcontroller',
            executable='controller_node',
            name='controller_node',
            output='screen',
            namespace='ns1_368540'
        ),
        Node(
            package='turtlesim',
            executable='turtlesim_node',
            name='sim2',
            output='screen',
            namespace='ns2_368540'
        ),
        Node(
            package='mybestcontroller',
            executable='controller_node',
            name='controller_node',
            output='screen',
            namespace='ns2_368540'
        )
    ])