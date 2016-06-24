# SolarSystemOrbiter
Plot the orbits of the planets in our Solar System and calculate the [Hohmann Transfer Orbits](https://en.wikipedia.org/wiki/Hohmann_transfer_orbit) to transfer your rocket ship from one to the other and back. Plot your route and create a travel movie to show your family and friends!

![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/htm/venus_mars.png?raw=true)

### How to install
Download the .zip or clone the repository, and you're done! 
The script is written in *Python 3.5* . It requires the *numpy*, *matplotlib*, *tkinter*, *imageio*, and the beautiful [*minibar*](https://github.com/canassa/minibar) package. All these packages can be acquired with the pip installer
`pip install *package*`

### How to use it

The script is run using
`python sso.py`

This interface will then pop up:

![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/htm/interface.png?raw=true)

You choose the planets that shall be plotted and then select the origin and desitnation of your travel.

Selecting the 'Acceleration at Apohelion' button triggers a second acceleration as soon as you have reachead your destination, meaning that the rocket ship stays in the new orbit. 

Clicking the 'Get Rocket Trajectory' button calculates the required time (number of integration steps) to reach your target, and the required planetary offset in degree so you don't miss Mars or the other planets once you reach the orbit. You can change the values and expirement what would happen, if..!

The 'Plot' button triggers the integration of the orbits of the planets and the transfer and opens a plot for you. The 'Make Movie' button saves 50 plots and creates a .gif in the folder that you specify in the entry field to the left.

Some status information is shown in the command line while you run the script!
![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/htm/progress.png?raw=true)
### What it does
The planet and transfer orbits are calculated using the [leap-frog integration scheme](https://en.wikipedia.org/wiki/Leapfrog_integration). Several assumptions are made, mostly in-plane, circular orbits of the planets and the rocket ship around the Sun. Only the gravitational field of the Sun is regarded.

Here is the math behind the calculation of the HTM and an examplary calculation for a Earth-Mars transfer:
![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/htm/maths.png?raw=true)


### To-Do
* Order of images in GIF appears to be messed up
* Include Planet IX 