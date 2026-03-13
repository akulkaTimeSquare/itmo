import glob
from setuptools import find_packages, setup

package_name = 'my_best_robot_simulation_setup'

setup(
    name=package_name,
    version='0.0.0',
    packages=find_packages(exclude=['test']),
    data_files=[
        ('share/ament_index/resource_index/packages',
            ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
        ('share/' + package_name + "/launch", ["launch/simulation_final.launch.py"]),
        ('share/' + package_name + "/urdf", glob.glob("urdf/*")),
        ('share/' + package_name + "/worlds", glob.glob("worlds/*")),
        ('share/' + package_name + "/config", glob.glob("config/*")),
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
            'controller = my_best_robot_simulation_setup.controller:main',
            'data_proc = my_best_robot_simulation_setup.data_proc:main',
        ],
    },
)
