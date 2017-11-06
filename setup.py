try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup



setup(
    name="SolarSystemOrbiter",
    url="https://github.com/madoee/SolarSystemOrbiter",
    author="Max Mahlke",
    version="2.0",
    packages=['SolarSystemOrbiter'],
    author_email='m_mahlke@yahoo.de',
    description="Plot orbits of planets and calculate Hohmann Tranfer Orbits",
    long_description=open("README.md").read(),
    install_requires=[
        "matplotlib",
        "numpy",
    ]
)
