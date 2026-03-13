from glob import glob
from setuptools import find_packages, setup

package_name = 'my_best_robot_simulation_control'

setup(
    name=package_name,
    version='0.0.0',
    packages=find_packages(exclude=['test']),
    data_files=[
        ('share/ament_index/resource_index/packages',
            ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
        ('share/' + package_name + "/launch", ["launch/simulation_control.launch.py"]),
        ('share/' + package_name + "/urdf", glob("urdf/*")),
        ('share/' + package_name + "/worlds", glob("worlds/*")),
        ('share/' + package_name + "/config", glob("config/*")),
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
            'controller = my_best_robot_simulation_control.controller:main',
        ],
    },
)
