import numpy as np


def hohmann(nsteps, d_o, d_d, second):
    # d_o is distance of origin to sun in AU
    # d_d is distance of destination to sun in AU
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

    AU = 150e9
    G = 6.67e-11				# gravitational constant
    M = 1.998e30					# mass of Sun

    angle = 0
    a = angle * np.pi / 180
    d_init = d_o * AU
    r_1 = d_d * AU
    v = np.sqrt(G*M/d_init) + np.sqrt(G*M/d_init) * (np.sqrt(2*r_1 / (r_1 + d_init)) - 1)

    x_0 = d_init*np.sin(a)
    y_0 = -d_init*np.cos(a)
    v_x_0 = v*np.cos(a)
    v_y_0 = v*np.sin(a)

    delta = 2 * 3.141597 * 1*AU / np.sqrt(G*M/AU) / 365 / 15     # Step size is 1/15 Earth day

    # Taylor approximation of velocities at n=1/2
    v_x_0 -= G*M*x_0/r(x_0, y_0)**3 * delta/2
    v_y_0 -= G*M*y_0/r(x_0, y_0)**3 * delta/2

    second_impulse = False

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
        if second:
            if y_0 >= r_1:
                if not second_impulse:
                    print('Impulse at Apohelion')
                    v_x_0 -= np.sqrt(G*M/d_d/AU) * (1 - np.sqrt(2*d_o*AU / (d_o*AU + d_d*AU)))
                    # v_y_0 += np.sqrt(G*M/d_d/AU) * (1 - np.sqrt(2*d_o*AU / (d_o*AU + d_d*AU)))
                    v_x_1 -= np.sqrt(G*M/d_d/AU) * (1 - np.sqrt(2*d_o*AU / (d_o*AU + d_d*AU)))
                    # v_y_1 += np.sqrt(G*M/d_d/AU) * (1 - np.sqrt(2*d_o*AU / (d_o*AU + d_d*AU)))
                    second_impulse = True
        yield x_0 / AU, y_0 / AU


def leap_frog(nsteps, distance, angle):
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

    AU = 150e9
    G = 6.67e-11				# gravitational constant
    M = 1.998e30					# mass of Sun

    # Use distance and anlge to calculate inital positions x_0, y_0 and velocities
    d = distance * AU
    a = angle * np.pi / 180
    v = np.sqrt(G*M/d)

    x_0 = d*np.sin(a)
    y_0 = -d*np.cos(a)

    v_x_0 = v*np.cos(a)
    v_y_0 = v*np.sin(a)

    delta = 2 * 3.141597 * 1*AU / np.sqrt(G*M/AU) / 365 / 15     # Step size is 1 Earth day
    # Taylor approximation of velocities at n=1/2
    v_x_0 = v_x_0 - G*M*x_0/r(x_0, y_0)**3 * delta/2
    v_y_0 = v_y_0 - G*M*y_0/r(x_0, y_0)**3 * delta/2

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
        yield x_0 / AU, y_0 / AU
