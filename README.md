# SolarSystemOrbiter
Plot the orbits of the planets in our Solar System and calculate the [Hohmann Transfer Orbits](https://en.wikipedia.org/wiki/Hohmann_transfer_orbit) and [Gravity Assist Maneuvers](https://en.wikipedia.org/wiki/Gravity_assist) to transfer your rocket ship from one to the other and back.

![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/earth_mars.png?raw=true)

## Updated Version - Now with Gravity Assist Maneuvers
The code is now much simpler to read and understand, using a "semi-class" structure (i.e. a dictionary) and generators for computation of the orbits. Now the spacecraft "knows" where the planets are at all times, allowing for Gravity Assist Maneuvers.

At the moment, there is no suggested "Planet Offset and Initial Impulse"-Calculation function to really use GAMs, this will be the next step. For now, you can use the *GAM* option using HTs to get your spacecraft close to the swing-by planet. The program will load a nice configuration at start-up. You can adjust the *Planet Offset* value to experiment with the GAM parameters. In general, the GAMs work better for larger HT distances, as the uncertanties lead to the spacecraft slightly missing the planet. For close trips, the calculations lead to the spacecraft "crashing into" the planet (the planets are point sources in this simulation though).

For GAMs, the accuracy has to be increased. This leads to longer integration times.

The instructions below are a little outdated, but summarize the main functionality. I will get back to updating this soon.

### How to install
Download the .zip or clone the repository, and you're done!
The script is written in Python 3, but also Python 2 compatible. It requires the *numpy*, *matplotlib*, *seaborn*, *tkinter*, and *imageio*. All these packages can be acquired with the pip installer

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
The orbits of the planets reflect their actual distances and eccentricities. If a Hohmann Transfer Orbit is calculated, the script sets the eccentricity of destination and origin orbit to 0, as the varying speed of the planets on elliptical orbits messes up the trajectory calculation (will be implemented later).

Here is the math behind the calculation of the HTM and an examplary calculation for a Earth-Mars transfer:
![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/maths.png?raw=true)

We disregard any perturbations due to third bodies, non-spherical inhomogenous planets, solar radiation, ...
In the end, all we do is accelerate our rocket ship at the start and see how its trajectory is influenced by the Sun's gravity field.

### To-Do
* Increase number of integration steps only for GAM calculation in SoI
* Include Planet IX
* Include *One-Tangent Burns*. Compare to HTs
* Combine HTs with phasing orbits
* Use eccentric anomaly to calculate HTs for elliptical orbits (fixed orientation between origin and target planet)
* Add real planet's positions for acutal dates. Look for possible next date for HT / GAM between planets.
* Long-term goal: Itineraries of Hohmann Transfers and GAMs. Comparison of energies and times required
