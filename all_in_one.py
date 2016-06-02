from __future__ import division
import numpy as np
import math

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
    a = angle * math.pi / 180
    d_init = d_o * AU
    r_1 = d_d * AU
    v = np.sqrt(G*M/d_init) + np.sqrt(G*M/d_init) * (np.sqrt(2*r_1 / (r_1 + d_init)) - 1)

    x_0 = d_init*np.sin(a)
    y_0 = -d_init*np.cos(a)
    v_x_0 = v*np.cos(a)
    v_y_0 = v*np.sin(a)



    #v_x_0 = 0
    #y_0 = 0
    # Initial velocity is orbital velocity of origin + impulse to change orbit. Impulse is negative or positive,
    # depending on relation between d_o and d_d
    #v_y_0 =     # 28 km / s



#    delta = 2 * 3.141597 * (x_0 + r_1 / 2) / v_y_0 / nsteps   # Now have an ellipse with perimeter x_0 + r_1 /2
    delta = 2 * 3.141597 * 1*AU / np.sqrt(G*M/AU) / 365 / 15     # Step size is 1 Earth day

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

    def r(x,y):
        # Calculate Earth-Sun distance
        r_0 = np.sqrt(x**2 + y**2)
        return r_0

    AU = 150e9
    G = 6.67e-11				# gravitational constant
    M = 1.998e30					# mass of Sun

    # Use distance and anlge to calculate inital positions x_0, y_0 and velocities
    d = distance * AU
    a = angle * math.pi / 180
    v = np.sqrt(G*M/d)

    x_0 = d*np.sin(a)
    y_0 = -d*np.cos(a)

    v_x_0 = v*np.cos(a)
    v_y_0 = v*np.sin(a)


    #x_0 = x_0 * AU
    #v_x_0 = 0
    #y_0 = 0
    #v_y_0 = np.sqrt(G*M/x_0)	# 28 km / s


    delta = 2 * 3.141597 * 1*AU / np.sqrt(G*M/AU) / 365      # Step size is 1 Earth day
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

#
# # Lists to store orbit coordinates
# planet_x = []
# planet_y = []
# hohmann_x = []
# hohmann_y = []
#
# x2 = []
# y2 = []
# # counter for images
# i = 0
# generator = hohmann(7*37, 1., 1.5)
# generator2 = leap_frog(7*37, 1.5, 44.)
# for x,y in leap_frog(7*37, 1., 0.):
#     planet_x.append(x)
#     planet_y.append(y)
#     a, b = next(generator)
#     c, d = next(generator2)
#     x2.append(c)
#     y2.append(d)
#     hohmann_x.append(a)
#     hohmann_y.append(b)
#     i += 1
#     if i % 10 == 0:
#         fig = plt.Figure()
#         plt.plot(planet_x, planet_y, 'b-')
#         plt.plot(x, y, 'o', ms=8)
#         plt.plot(x2, y2, 'b-')
#         plt.plot(c, d, 'o', ms=8)
#         # Sun
#         plt.plot(0, 0, 'yo', ms=10)
#         plt.plot(hohmann_x, hohmann_y, 'r-')
#         plt.plot(a,b, 'x', ms=5)
#         # Set plot limits to largest radius + x
#         plt.grid()
#         plt.xlim(-1.7, 1.7)
#         plt.ylim(-1.7, 1.7)
#         plt.xlabel("x / AU")
#         plt.ylabel("y / AU")
#         plt.savefig('/Users/Max/Desktop/CA/' + str(i) + '.png')
#         fig.clear()
#         plt.cla()
#         plt.clf()
#
# images = []
# image_files = ['/Users/Max/Desktop/CA/' + f for f in os.listdir('/Users/Max/Desktop/CA/' + '/') if os.path.isfile(os.path.join('/Users/Max/Desktop/CA/', f)) and f[-4:] == '.png']
# for filenam in image_files:
#     images.append(imageio.imread(filenam))
#     imageio.mimsave('/Users/Max/Desktop/epic.gif', images)        # Save object dictionary entry to text file
#
# from __future__ import division
# import numpy as np
# import math
# import matplotlib.pyplot as plt
#
# def hohmann(nsteps, d_o, d_d):
#     # d_o is distance of origin to sun in AU
#     # d_d is distance of destination to sun in AU
#     def step_x():
#         # Calculate step in x, return incremented x
#         x_1 = x_0 + delta * v_x_0
#         return x_1
#
#     def step_y():
#         # Calculate step in y, return incremented y
#         y_1 = y_0 + delta * v_y_0
#         return y_1
#
#     def step_v_x():
#         # Calculate step in v_x, return incremented v_x
#         v_x = v_x_0 - delta * G*M*x_1/r(x_1, y_1)**3
#         return v_x
#
#
#     def step_v_y():
#         # Calculate step in v_y, return incremented v_y
#         v_y = v_y_0 - delta * G*M*y_1/r(x_1, y_1)**3
#         return v_y
#
#     def r(x, y):
#         # Calculate Earth-Sun distance
#         r_0 = np.sqrt(x**2 + y**2)
#         return r_0
#
#     AU = 150e9
#     G = 6.67e-11				# gravitational constant
#     M = 1.998e30					# mass of Sun
#
#     angle = 0
#     a = angle * math.pi / 180
#     d_init = d_o * AU
#     r_1 = d_d * AU
#     v = np.sqrt(G*M/d_init) + np.sqrt(G*M/d_init) * (np.sqrt(2*r_1 / (r_1 + d_init)) - 1)
#
#     x_0 = d_init*np.sin(a)
#     y_0 = -d_init*np.cos(a)
#     v_x_0 = v*np.cos(a)
#     v_y_0 = v*np.sin(a)
#
#
#
#     #v_x_0 = 0
#     #y_0 = 0
#     # Initial velocity is orbital velocity of origin + impulse to change orbit. Impulse is negative or positive,
#     # depending on relation between d_o and d_d
#     #v_y_0 =     # 28 km / s
#
#
#
# #    delta = 2 * 3.141597 * (x_0 + r_1 / 2) / v_y_0 / nsteps   # Now have an ellipse with perimeter x_0 + r_1 /2
#     delta = 2 * 3.141597 * 1*AU / np.sqrt(G*M/AU) / 365      # Step size is 1 Earth day
#
#     # Taylor approximation of velocities at n=1/2
#     v_x_0 -= G*M*x_0/r(x_0, y_0)**3 * delta/2
#     v_y_0 -= G*M*y_0/r(x_0, y_0)**3 * delta/2
#
#     second_impulse = False
#
#     x = []
#     y = []
#
#     # Integration
#     for step in range(0, nsteps):
#         x.append(x_0 / AU)
#         y.append(y_0 / AU)
#         x_1 = step_x()
#         y_1 = step_y()
#         v_x_1 = step_v_x()
#         v_y_1 = step_v_y()
#         x_0 = x_1
#         y_0 = y_1
#         v_x_0 = v_x_1
#         v_y_0 = v_y_1
#         if step % 10 == 0:
#             plt.plot(x, y)
#             plt.savefig('/Users/Max/Desktop/CA/' + str(step) + '.png')
#         # Check if we have reached Target yet
#         #if y_0 >= r_1:
#         #    if not second_impulse:
#         #        v_y_0 += np.sqrt(G*M/d_o/AU) * (1 - np.sqrt(2*d_d*AU / (d_o*AU + d_d*AU)))
#         #        v_y_1 += np.sqrt(G*M/d_o/AU) * (1 - np.sqrt(2*d_d*AU / (d_o*AU + d_d*AU)))
#         #        second_impulse = True
#
#     return x, y
