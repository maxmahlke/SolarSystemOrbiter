# SolarSystemOrbiter - 2.0
Plot the orbits of the planets in our Solar System and calculate the [Hohmann Transfer Orbits](https://en.wikipedia.org/wiki/Hohmann_transfer_orbit) to transfer your rocket ship from one planet to the other and back.

![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/progress.png?raw=true)

### How to install
Download the .zip or clone the repository, and you're done!
The script is written in Python 3, but also Python 2 compatible. It requires the *numpy*, *matplotlib*, and *tkinter* packages. All these packages can be acquired with the pip installer

`pip install *package*`

or by executing the install script

`[sudo] python3 setup.py install`

### How to use it

The script is run using
`python sso.py`

This interface will then pop up:

![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/interface.png?raw=true)

On the top you choose the planets you want to include in the simulation. If you want to calculate a Hohmann Transfer Orbit, choose the origin and destination in the lists below.

For the transfer orbit to be plotted, you have to select the *Plot Hohmann Transfer Orbit* checkbox. Selecting the *Orbit Insertion* checkbox triggers a second acceleration as soon as you have reached your destination (defined by the sphere of influence of the destination planet). The spacecraft now follows the target planet's orbit.

The *Timestep* and *Integration Time* textboxes set the simulation parameters. A finer (smaller) timestep will lead to more accurate results, while increasing the computation time. The *Suggest Simulation Parameters* button proposes values based on the planets included in the simulation and the parameters of the Hohmann Transfer orbit.

The *Plot* button triggers the integration of the orbits of the planets and the transfer and opens a plot for you. The *Animate* button shows you the results of the simulation as they are calculated, with additional information superimposed. This animation can be saved to the local directory using the *Save as mp4* button. A file called *SolarSystemSimulation.mp4* is created.


### What it does
The planet and transfer orbits are calculated using the [leap-frog integration scheme](https://en.wikipedia.org/wiki/Leapfrog_integration). Several assumptions are made: in-plane orbits of the planets and the rocket ship around the Sun, only the gravitational field of the Sun is regarded.
The orbits of the planets reflect their actual distances and eccentricities. If a Hohmann Transfer Orbit is calculated, the script sets the eccentricity of destination and origin orbit to 0, as the varying speed of the planets on elliptical orbits messes up the trajectory calculation (will be implemented later).

Here is the math behind the calculation of the HTM and an exemplary calculation for a Earth-Mars transfer:
![alt-tag](https://github.com/madoee/SolarSystemOrbiter/blob/master/SolarSystemOrbiter/htm/maths.png?raw=true)

We disregard any perturbations due to third bodies, non-spherical inhomogeneous planets, solar radiation, ...
In the end, all we do is accelerate our rocket ship at the start and see how its trajectory is influenced by the Sun's gravity field.

### To-Do
* Increase number of integration steps only for GAM calculation in SoI
* Include Planet IX
* Include *One-Tangent Burns*. Compare to HTs
* Combine HTs with phasing orbits
* Use eccentric anomaly to calculate HTs for elliptical orbits (fixed orientation between origin and target planet)
* Add real planet's positions for actual dates. Look for possible next date for HT / GAM between planets.
* Long-term goal: Itineraries of Hohmann Transfers and GAMs. Comparison of energies and times required
