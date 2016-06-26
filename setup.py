try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup



setup(
    name="SolarSystemOrbiter",
    url="https://github.com/madoee/SolarSystemOrbiter",
    author="Max @madoee",
    description="Plot orbits of our planets",
    long_description=open("README.md").read(),
    install_requires=[
        "matplotlib",
        "numpy",
        "imageio",
        "minibar",
    ]
)
