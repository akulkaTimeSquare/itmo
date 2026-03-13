from setuptools import find_packages, setup

package_name = 'my_best_model'

setup(
    name=package_name,
    version='0.0.0',
    packages=find_packages(exclude=['test']),
    data_files=[
        ('share/ament_index/resource_index/packages',
            ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
        ('share/' + package_name + "/launch", ["launch/rviz.launch.py"]),
        ('share/' + package_name + "/config", ["config/rviz.rviz"]),
        ('share/' + package_name + "/urdf", ["urdf/robot_model.xacro",
                                              "urdf/materials.xacro",
                                              "urdf/wheel.xacro"]),
    ],
    install_requires=['setuptools'],
    zip_safe=True,
    maintainer='yorides',
    maintainer_email='akulkaTimeSquare@gmail.com',
    description='TODO: Package description',
    license='Apache-2.0',
    extras_require={
        'test': [
            'pytest',
        ],
    },
    entry_points={
        'console_scripts': [
        ],
    },
)
