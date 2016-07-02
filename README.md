# SolarSystemOrbiter
Plot the orbits of the planets in our Solar System and calculate the [Hohmann Transfer Orbits](https://en.wikipedia.org/wiki/Hohmann_transfer_orbit) to transfer your rocket ship from one to the other and back.

![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/earth_mars.png?raw=true)

### How to install
Download the .zip or clone the repository, and you're done!
The script is written in Python 3, but also Python 2 compatible. It requires the *numpy*, *matplotlib*, *seaborn*, *tkinter*, *imageio*, and the [*minibar*](https://github.com/canassa/minibar) package. All these packages can be acquired with the pip installer

`pip install *package*`  

or by executing the install script  

`[sudo] python3 setup.py install`

### How to use it

The script is run using
`python sso.py`

This interface will then pop up:

![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/interface.png?raw=true)

On the top you choose the planets you want to include in the simulation. If you want to calculate a Hohmann Transfer Orbit, choose the origin and destination in the lists below.

For the transfer orbit to be plotted, you have to select the 'Plot Hohmann Transfer Orbit' checkbox. Selecting the 'Acceleration at Apohelion' button triggers a second acceleration as soon as you have reachead your destination, meaning that the rocket ship stays in the new orbit.

Clicking the 'Get Rocket Trajectory' button calculates the required time (number of integration steps) to reach your target, and the required planetary offset in degree so you don't miss Mars or the other planets once you reach the orbit. You can change the values and expirement what would happen, if..!

The 'Plot' button triggers the integration of the orbits of the planets and the transfer and opens a plot for you. The 'Make Movie' button saves 50 plots and creates a .gif in the folder that you specify in the entry field to the left.

Some status information is shown in the command line while you run the script!
![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/progress.png?raw=true)


### What it does
The planet and transfer orbits are calculated using the [leap-frog integration scheme](https://en.wikipedia.org/wiki/Leapfrog_integration). Several assumptions are made: in-plane orbits of the planets and the rocket ship around the Sun, only the gravitational field of the Sun is regarded.
The orbits of the planets reflect their actual distances and eccentricities. If a Hohmann Transfer Orbit is calculated, the script sets the eccentricity of destination and origin orbit to 0, as the varying speed of the planets on elliptical orbits messes up the trajectory calculation.

Here is the math behind the calculation of the HTM and an examplary calculation for a Earth-Mars transfer:
![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/maths.png?raw=true)


### To-Do
* Add names of planet to orbits
* Order of images in GIF appears to be messed up
* Include Planet IX
* Include Gravity Assist Orbits
