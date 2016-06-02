import numpy as np


def suggest(d_o, d_d, T_d):
    # Calculates travel time of rocket and required angular offset of planets

    # Travel time
    # Using Keplers third law units of 1 AU and 1 Earth year: T^2 = a^3

    a = (d_o + d_d) / 2
    T_squared = a**3
    T = np.sqrt(T_squared)
    # Only half-way required
    T /= 2

    # Calculate required angular offset
    off = T/T_d * 360
    # Return T in Earth Weeks
    return T*52, 180-off


