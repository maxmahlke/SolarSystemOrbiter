import numpy as np


def leap_frog(nsteps, semi_major, eccentricity, angle):
    def step_x():
        # Calculate step in x, return incremented x
        x_1 = x_0 + delta * v_x_0
        return x_1

    def step_y():
        # Calculate step in y, return incremented y
        y_1 = y_0 + delta * v_y_0
        return y_1

    def step_v_x():
        # Calculate step in v_x, return incremented v_x
        v_x = v_x_0 - delta * G*M*x_1/r(x_1, y_1)**3
        return v_x

    def step_v_y():
        # Calculate step in v_y, return incremented v_y
        v_y = v_y_0 - delta * G*M*y_1/r(x_1, y_1)**3
        return v_y

    def r(x, y):
        # Calculate Earth-Sun distance
        r_0 = np.sqrt(x**2 + y**2)
        return r_0

    AU = 149.6e9
    G = 6.67e-11				# gravitational constant
    M = 1.998e30					# mass of Sun

    # Use distance and anlge to calculate inital positions x_0, y_0 and velocities
    distance = semi_major * (1. - eccentricity)
    d = distance * AU
    theta = angle * np.pi / 180.
    v = np.sqrt(G*M*(2./d - 1./(AU*semi_major)))

    x_0 = d*np.sin(theta)
    y_0 = -d*np.cos(theta)

    v_x_0 = v*np.cos(theta)
    v_y_0 = v*np.sin(theta)

    delta = 2. * 3.141597 * AU / np.sqrt(G*M/AU) / 365. / 15.      # Step size is 1 / 15 Earth day
    # Taylor approximation of velocities at n=1/2
    v_x_0 = v_x_0 - G*M*x_0/r(x_0, y_0)**3 * delta/2.
    v_y_0 = v_y_0 - G*M*y_0/r(x_0, y_0)**3 * delta/2.

    x = []
    y = []

    # Integration
    for step in range(0, nsteps):
        x.append(x_0 / AU)
        y.append(y_0 / AU)
        x_1 = step_x()
        y_1 = step_y()
        v_x_1 = step_v_x()
        v_y_1 = step_v_y()
        x_0 = x_1
        y_0 = y_1
        v_x_0 = v_x_1
        v_y_0 = v_y_1

    return x, y
