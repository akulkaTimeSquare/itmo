import os
from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.actions import IncludeLaunchDescription
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import LaunchConfiguration, PathJoinSubstitution, Command
from launch_ros.actions import Node
from launch_ros.substitutions import FindPackageShare


def generate_launch_description():
    pkg_ros_gz_sim = get_package_share_directory('ros_gz_sim')
    pkg_lab = get_package_share_directory('my_best_robot_simulation_setup')
    gui_config_path = os.path.join(pkg_lab, 'config', 'gaz.config')
    rviz_config_path = os.path.join(pkg_lab, 'config', 'rviz.rviz')
    
    xacro_file = os.path.join(pkg_lab, 'urdf', 'robot_model.xacro')
    world_path = PathJoinSubstitution([
        FindPackageShare('my_best_robot_simulation_setup'),
        'worlds','demo.world'
    ])

    # Gazebo
    gz_sim = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(pkg_ros_gz_sim, 'launch', 'gz_sim.launch.py')),
        launch_arguments={
            'gz_args': [world_path, ' -r', ' --gui-config ', gui_config_path],
        }.items()
    )

    # Rviz
    rviz = Node(
        package='rviz2',
        executable='rviz2',
        parameters=[{'use_sim_time': True}],
        arguments=['-d', rviz_config_path],
        output='screen',
    )

    # Robot State Publisher
    robot_state_publisher = Node(
        package='robot_state_publisher',
        executable='robot_state_publisher',
        name='robot_state_publisher',
        output='both',
        parameters=[
            {'use_sim_time': True},
            {'robot_description': Command(['xacro ', LaunchConfiguration('urdf_model')])},
        ],
    )

    # Spawn robot
    spawn = Node(
        package='ros_gz_sim',
        executable='create',
        arguments=[
            '-name', 'diff_drive',
            '-topic', 'robot_description',
            '-x', '0', '-y', '0', '-z', '0.2'
        ],
        output='screen',
    )

    # Bridge for ROS-GZ communication
    bridge = Node(
        package='ros_gz_bridge',
        executable='parameter_bridge',
        parameters=[{
            'config_file': os.path.join(pkg_lab, 'config', 'ros_gz_bridge.yaml'),
            'qos_overrides./tf_static.publisher.durability': 'transient_local'
        }],
        output='screen',
    )

    # Controller for /cmd bridge
    spawn_controller = Node(
        package="controller_manager",
        executable="spawner",
        name=f"spawn_controller_wheel_controller",
        arguments=["velocity_controller"],
        output="screen",
    )

    # Controller node
    cmd_vel_controller = Node(
        package='my_best_robot_simulation_setup',
        executable='controller',
        name='controller',
        output='screen',
    )

    # Data processing node
    data_proc_node = Node(
        package='my_best_robot_simulation_setup',
        executable='data_proc',
        name='data_proc',
        output='screen',
    )

    return LaunchDescription([
        gz_sim,
        DeclareLaunchArgument(
            'urdf_model',
            default_value=xacro_file,
            description='Full path to robot the Xacro file'
        ),
        robot_state_publisher,
        rviz,
        spawn,
        bridge,
        spawn_controller,
        cmd_vel_controller,
        data_proc_node,
    ])