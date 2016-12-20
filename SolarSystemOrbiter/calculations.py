import numpy as np
import minibar
# In case we want to plot the GAM trajectory
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
plt.style.use('seaborn-dark')

def soi(d, m):
    # Calculates Sphere of Influcene according to
    # r = d * (m/M)^(2/5) http://www.bogan.ca/orbits/gravasst/gravasst.html
    M = 1.998e30	               # mass of Sun
    return d*np.power(m/M, 0.4)     # in AU

def step_position(current_position, current_velocity, delta):
	# Calculate change in x or y coordinate, return updated coordinate
    updated_position = current_position + delta * current_velocity
    return updated_position

def r(x, y):
    # Calculate spacecraft - central-body distance
    r_0 = np.sqrt(x**2 + y**2)
    return r_0

def step_velocity(current_velocity, current_position_parallel, current_position_orthogonal, delta, mass_of_central_body):
    # Calculate step in velocity component, return updated velocity
    # Component depends on corresponding position component (here called parallel and orthogonal) and distance to central body
    # Note the order of the passed function variables is used to differentiate v_x and v_y
    G = 6.67e-11                # gravitational constant
    updated_velocity = current_velocity - delta * G*mass_of_central_body*current_position_parallel/r(current_position_parallel, current_position_orthogonal)**3
    return updated_velocity

def hohmann(steps, origin, target, mode, movie):
    # mode can be either HT with second impulse (orbit injection) or GAM or none
 
    print(u'\nLaunching Spacecraft..\n')

    G = 6.67e-11                # gravitational constant
    M = 1.998e30                # mass of Sun

    angle = 0.
    a = angle * np.pi / 180
    
    # d_o is distance of origin to sun in AU
    # d_d is distance of destination to sun in AU
    d_o = origin['semi_major'] * (1 - origin['eccentricity'])  # Leave at perihelion
    d_d = target['semi_major']  # Arrive at apohelion
    d_init = d_o
    r_1 = d_d

    v = np.sqrt(G*M/d_init) + np.sqrt(G*M/d_init) * (np.sqrt(2*r_1 / (r_1 + d_init)) - 1)

    x_0 = d_init*np.sin(a)
    y_0 = -d_init*np.cos(a)
    v_x_0 = v*np.cos(a)
    v_y_0 = v*np.sin(a)

    delta = 2 * 3.141597 /4 / np.sqrt(G*M) / 3650 / 15     # Step size is 1 Earth day

    # Taylor approximation of velocities at n=1/2
    v_x_0 -= G*M*x_0/r(x_0, y_0)**3 * delta/2
    v_y_0 -= G*M*y_0/r(x_0, y_0)**3 * delta/2


    # Maneuvers at target planet. We can do nothing, perform an GAM or use a second impulse
    # to insert the spacecraft into the orbit of the target planet. If GAM, the SoI is relevant.
    sphere_of_influence = r_1*0.0001

    if mode == 2: # GAM
       	# Caluclate SoI to define points where GAM begins / ends
        sphere_of_influence = soi(r_1, target['mass'])
        # variable for the instantaneous impulse when in SoI of swing-by planet
        gam_impulse = False

    elif mode == 1: # Second Impulse
    	second_impulse = False

    x = []
    y = []

    # Integration
    #for step in range(0, steps):
    while True:
        x.append(x_0)
        y.append(y_0)
        x_1 = step_position(x_0, v_x_0, delta)
        y_1 = step_position(y_0, v_y_0, delta)
        v_x_1 = step_velocity(v_x_0, x_1, y_1, delta, M)
        v_y_1 = step_velocity(v_y_0, y_1, x_1, delta, M)
        x_0 = x_1
        y_0 = y_1
        v_x_0 = v_x_1
        v_y_0 = v_y_1
        target_x, target_y = yield

        if mode != 0:  # If we perform GAM or Orbit injection
            # See if space-craft is close to target planet 
            # close being sphere of influcence for GAMs
            distance_to_target = np.sqrt((x_0-target_x)**2+(y_0-target_y)**2)
            if distance_to_target <= sphere_of_influence:
                # We have entered the SoI. We can now either perform a GAM or insert our spacecraft
                # into the planet's orbit. Both will result in an acceleration of the spacecraft.
                # If we do nothing, the spacecraft will simply follow its current orbit (we disregard that it entered the SoI)
                if mode == 2: # GAM
                    # If GAM has not been performed yet
                    if not gam_impulse:
                        print(u'\U0001f680 : Reached Sphere of Influence of Target Planet.\n    Performing Gravity Assist Maneuver..')
                        v_x_0, v_y_0 = GAM((x_1, y_1), (target_x, target_y), (v_x_1, v_y_1), target, movie)
                        #plt.plot((x_1, target_x), (y_1, target_y))
                        gam_impulse = True # impulse has been provided now
                elif mode == 1: # Orbit injection
                    # If we want to insert the spacecraft into the orbit of the target planet
                    if not second_impulse:
                        # If insertion has not been performed yet
                        print(u'\U0001f680 : Reached Sphere of Influence of Target Planet.\n    Performing Orbit Insertion Maneuver..')
                        v_x_0 -= np.sqrt(G*M/d_d) * (1 - np.sqrt(2*d_o / (d_o + d_d)))
                        v_x_1 -= np.sqrt(G*M/d_d) * (1 - np.sqrt(2*d_o / (d_o + d_d)))
                        second_impulse = True
        yield (x_0, y_0)


def planets(steps, semi_major, eccentricity, angle):

    G = 6.67e-11    			# gravitational constant
    M = 1.998e30					# mass of Sun

    # Use distance and anlge to calculate inital positions x_0, y_0 and velocities
    distance = semi_major * (1. - eccentricity)
    d = distance
    theta = angle * np.pi / 180.
    v = np.sqrt(G*M*(2./d - 1./semi_major))

    x_0 = d*np.sin(theta)
    y_0 = -d*np.cos(theta)

    v_x_0 = v*np.cos(theta)
    v_y_0 = v*np.sin(theta)

    delta = 2 * 3.141597 /4 / np.sqrt(G*M) / 3650 / 15     # Step size is 1 Earth day
    # Taylor approximation of velocities at n=1/2
    v_x_0 = v_x_0 - G*M*x_0/r(x_0, y_0)**3 * delta/2
    v_y_0 = v_y_0 - G*M*y_0/r(x_0, y_0)**3 * delta/2

    x = []
    y = []

    # Integration
    for step in range(0, steps):
        x.append(x_0)
        y.append(y_0)
        x_1 = step_position(x_0, v_x_0, delta)
        y_1 = step_position(y_0, v_y_0, delta)
        v_x_1 = step_velocity(v_x_0, x_1, y_1, delta, M)
        v_y_1 = step_velocity(v_y_0, y_1, x_1, delta, M)

        x_0 = x_1
        y_0 = y_1
        v_x_0 = v_x_1
        v_y_0 = v_y_1
        yield x_0, y_0


def GAM(position_spacecraft, position_planet, velocity_spacecraft, target, movie):
	# We should now be either before or behind the swing-by planet.
	# We switch from heliocentric reference system to system of planet

	G = 6.67e-11    			# gravitational constant
	M = 1.998e30				# mass of Sun

	# x and y of spacecraft, heliocentric
	x_sc = position_spacecraft[0]
	y_sc = position_spacecraft[1]

	# x and y of planet, heliocentric
	x_p = position_planet[0]
	y_p = position_planet[1]

	# vx and vy of spacecraft, heliocentric
	vx_sc = velocity_spacecraft[0]	# Currently using HTs which lead us directly to the planet
	vy_sc = velocity_spacecraft[1]

	# vx and vy of planet, heliocentric
	v_p = np.sqrt(G*M/target['semi_major'])
	vx_p = -np.cos(np.arctan(x_p/y_p))*v_p
	vy_p = np.sin(np.arctan(x_p/y_p))*v_p

	# x and y of spacecraft, planet-centric
	x_sc = x_sc - x_p
	y_sc = y_sc - y_p

	# vx and vy of spacecraft, planet-centric
	vx_sc = vx_sc - vx_p
	vy_sc = vy_sc - vy_p

	# "Switch gravity" from Sun to swing-by planet.
	# Integration
	x_gam = []
	y_gam = []
	delta = 2 * 3.141597  / np.sqrt(G*M) / 3650 / 15
	M = 1.899e27 # mass of jupiter

	# calculate SoI
	sphere_of_influence  = soi(target['semi_major'], target['mass'])
	# calculate trajectory inside SoI
	if True:
		for step in range(0, int(1e5)):
			x_gam.append(x_sc)
			y_gam.append(y_sc)
			x_sc_1 = step_position(x_sc, vx_sc, delta)
			y_sc_1 = step_position(y_sc, vy_sc, delta)
			v_x_1 = step_velocity(vx_sc, x_sc_1, y_sc_1, delta, M)
			v_y_1 = step_velocity(vy_sc, y_sc_1, x_sc_1, delta, M)
			x_sc = x_sc_1
			y_sc = y_sc_1
			vx_sc = v_x_1
			vy_sc = v_y_1
			if x_sc**2+y_sc**2 > sphere_of_influence:
				# if we have left SoI
				break
		f = plt.figure(1)
		plt.plot(x_gam, y_gam, linestyle='--', color='Maroon')
		plt.plot(0., 0., marker='o', ms=8, color=target['color'])
		plt.xlabel('x / AU')
		plt.ylabel('y / AU')
		plt.grid()
		plt.title('Gravity Assist Maneuver')
		plt.axis('equal')
		if not movie:
			f.show()
		else:
			plt.savefig('gam.png')
			plt.gcf().clf()

	# Transform velocities back to heliocentric frame of reference
	# Planet's movement not considered...

	# x and y of spacecraft, helio-centric again
	x_sc += x_p
	y_sc += y_p

	# vx and vy of spacecraft, helio-centric again
	vx_sc += vx_p
	vy_sc += vy_p
	return vx_sc, vy_sc
	