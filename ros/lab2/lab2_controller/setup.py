from setuptools import find_packages, setup

package_name = 'lab2_controller'

setup(
    name=package_name,
    version='0.0.1',
    packages=find_packages(exclude=['test']),
    data_files=[
        ('share/ament_index/resource_index/packages',
            ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
        ('share/' + package_name + '/launch', 
            ['launch/service_server.launch.py','launch/service_client.launch.py',
             'launch/action_server.launch.py','launch/action_client.launch.py']),
    ],
    install_requires=['setuptools'],
    zip_safe=True,
    maintainer='yorides',
    maintainer_email='akulkaTimeSquare@gmail.com',
    description='TODO: Package description',
    license='Apache License 2.0',
    extras_require={
        'test': [
            'pytest',
        ],
    },
    entry_points={
        'console_scripts': [
            "service_server = lab2_controller.service_server:main",
            "service_client = lab2_controller.service_client:main",
            "action_server = lab2_controller.action_server:main",
            "action_client = lab2_controller.action_client:main",
        ],
    },
)
