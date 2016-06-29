#! /usr/bin/env python
try:
    import tkinter as tk
    from tkinter import ttk
    from _tkinter import TclError
except ImportError:  # Python 2
    import Tkinter as tk
    from _tkinter import TclError
    import ttk
from calc import all_in_one
from calc import hohmann as hm
from calc import leap_frog as lf
import imageio
import os
import numpy as np
import minibar
import matplotlib.pyplot as plt
plt.style.use('seaborn-dark')


class App:
    def __init__(self, master):
        self.master = master

        # Read-in catalogs architecture
        self.planet_frame = ttk.LabelFrame(master, text='Choose Planets to Plot')
        self.planet_frame.grid(row=0, columnspan=5, sticky='EW')

        # Planets
        self.mercury = tk.BooleanVar()
        self.venus = tk.BooleanVar()
        self.earth = tk.BooleanVar()
        self.mars = tk.BooleanVar()
        self.jupiter = tk.BooleanVar()
        self.saturn = tk.BooleanVar()
        self.uranus = tk.BooleanVar()
        self.neptune = tk.BooleanVar()

        self.plot_mercury = ttk.Checkbutton(self.planet_frame, text='Mercury', variable=self.mercury)
        self.plot_mercury.grid(row=0, column=0, sticky='EW', padx=5, pady=5)
        self.plot_venus = ttk.Checkbutton(self.planet_frame, text='Venus', variable=self.venus)
        self.plot_venus.grid(row=1, column=0, sticky='EW', padx=5, pady=5)
        self.plot_earth = ttk.Checkbutton(self.planet_frame, text='Earth', variable=self.earth)
        self.plot_earth.grid(row=2, column=0, sticky='EW', padx=5, pady=5)
        self.plot_mars = ttk.Checkbutton(self.planet_frame, text='Mars', variable=self.mars)
        self.plot_mars.grid(row=3, column=0, sticky='EW', padx=5, pady=5)
        self.plot_jupiter = ttk.Checkbutton(self.planet_frame, text='Jupiter', variable=self.jupiter)
        self.plot_jupiter.grid(row=0, column=1, sticky='EW', padx=5, pady=5)
        self.plot_saturn = ttk.Checkbutton(self.planet_frame, text='Saturn', variable=self.saturn)
        self.plot_saturn.grid(row=1, column=1, sticky='EW', padx=5, pady=5)
        self.plot_uranus = ttk.Checkbutton(self.planet_frame, text='Uranus', variable=self.uranus)
        self.plot_uranus.grid(row=2, column=1, sticky='EW', padx=5, pady=5)
        self.plot_neptune = ttk.Checkbutton(self.planet_frame, text='Neptune', variable=self.neptune)
        self.plot_neptune.grid(row=3, column=1, sticky='EW', padx=5, pady=5)

        self.transfer_frame = ttk.LabelFrame(master, text='Choose Origin and Destination')
        self.transfer_frame.grid(row=4, columnspan=5, sticky='EW')

        self.plot_hohmann = tk.BooleanVar()

        self.plot_hohmann_check = ttk.Checkbutton(self.transfer_frame, text='Plot Hohmann Transfer Orbit', variable=self.plot_hohmann)
        self.plot_hohmann_check.grid(row=4, column=0, sticky='EW', padx=5, pady=5)

        self.origin = tk.Listbox(self.transfer_frame, exportselection=0)
        self.origin.grid(row=5, column=0, sticky='EW')
        planets = ['Neptune', 'Uranus', 'Saturn', 'Jupiter', 'Mars', 'Earth', 'Venus', 'Mercury']
        for planet in planets:
            self.origin.insert(0, planet)

        self.destination = tk.Listbox(self.transfer_frame, exportselection=0)
        self.destination.grid(row=5, column=2, sticky='EW')
        for planet in planets:
            self.destination.insert(0, planet)

        # Planet offset in degree
        self.offs = tk.DoubleVar()
        # Boolean for acceleration at apohelion
        self.second = tk.BooleanVar()
        # File path for output images
        self.save_p = tk.StringVar()
        self.save_p.set(os.getcwd() + '/MyBigTrip/')

        tk.Entry(self.transfer_frame, textvariable=self.offs).grid(row=7, column=0, sticky='EW')
        tk.Label(self.transfer_frame, text='Planet Offset in Degree').grid(row=7, column=2, sticky='E')
        self.duration = tk.DoubleVar()
        tk.Entry(self.transfer_frame, textvariable=self.duration).grid(row=8, column=0, sticky='EW')
        tk.Label(self.transfer_frame, text='Flight Duration in Earth Weeks').grid(row=8, column=2, sticky='E')

        tk.Checkbutton(self.transfer_frame, text='Acceleration at Apohelion', variable=self.second).grid(row=6, column=0, sticky='EW')

        self.calc = ttk.Button(master, text='Get Rocket Trajectory', command=self.calc)
        self.calc.grid(column=0, row=9, sticky='EW')

        self.plot = ttk.Button(master, text='Plot', command=self.plot)
        self.plot.grid(column=2, row=9, sticky='EW')

        tk.Entry(self.master, textvariable=self.save_p).grid(row=10, column=0)
        self.movie = ttk.Button(master, text='Make Movie', command=self.movie)
        self.movie.grid(column=2, row=10, sticky='EW')

    def calc(self):
        # Have to select origin and destination of Hohmann Transfer. If not,
        # this function just returns a message and does nothing
        try:
            a, b = self.origin.get(self.origin.curselection()), self.destination.curselection()
        except TclError:
            print('You have to select origin and destiantion for the Hohmann Transfer Orbit')
            return 0

        print('Calculating Rocket trajectory..')

        planets = {'Mercury': [self.mercury.get(), 0.387, 0.241, 0.], 'Venus': [self.venus.get(), 0.723, 0.615, 0.],
                   'Earth': [self.earth.get(), 1., 1., 0.], 'Mars': [self.mars.get(), 1.524,  1.88, 0.],
                   'Jupiter': [self.jupiter.get(), 5.203, 11.9, 0.], 'Saturn': [self.saturn.get(), 9.58, 29.5, 0.],
                   'Uranus': [self.uranus.get(), 19.20, 84, 0.], 'Neptune': [self.neptune.get(), 30.06, 164.79, 0.]}

        # Calculates travel time of rocket and required angular offset of planets
        d_o = planets[self.origin.get(self.origin.curselection())][1]
        d_d = planets[self.destination.get(self.destination.curselection())][1]
        T_d = planets[self.destination.get(self.destination.curselection())][2]
        # Travel time
        # Using Keplers third law in units of 1 AU and 1 Earth year: T^2 = a^3
        a = (d_o + d_d) / 2
        T_squared = a**3
        T = np.sqrt(T_squared)
        # Only half-way required
        T /= 2
        # Calculate required angular offset
        off = T/T_d * 360
        # Return T in Earth Weeks
        t = T * 52
        off = 180-off

        self.duration.set(round(t, 2))
        self.offs.set(round(off, 2))

    def plot(self):
        # Have to select origin and destination of Hohmann Transfer if hohmann should be plotted. If not,
        # this function just returns a message and does nothing
        if self.plot_hohmann.get():
            try:
                a, b = self.origin.get(self.origin.curselection()), self.destination.curselection()
            except TclError:
                print('You have to select origin and destiantion for the Hohmann Transfer Orbit')
                return 0

        # Planet: Plotting Boolean, semi-major axis in AU, eccentricity, orbital period, angle, color for plot
        planets = {'Mercury': [self.mercury.get(), 0.387, 0.205, 0.241, 0., 'Gold'], 'Venus': [self.venus.get(), 0.723, 0.007, 0.615, 0., 'Coral'],
                   'Earth': [self.earth.get(), 1., 0.017, 1., 0., 'DarkBlue'], 'Mars': [self.mars.get(), 1.524, 0.094, 1.88, 0., 'Crimson'],
                   'Jupiter': [self.jupiter.get(), 5.203, 0.049, 11.9, 0., 'orange'], 'Saturn': [self.saturn.get(), 9.58, 0.057, 29.5, 0., 'Khaki'],
                   'Uranus': [self.uranus.get(), 19.20, 0.046, 84, 0., 'Turquoise'], 'Neptune': [self.neptune.get(), 30.06, 0.011, 164.79, 0., 'RoyalBlue']}

        nsteps = 7 * 15 * self.duration.get()
        if self.plot_hohmann.get():
            planets[self.destination.get(self.destination.curselection())][4] = self.offs.get()

        # Set planet offsets to random values, unless desitnation or target planet
        # Hohmann orbit optional
        for planet in planets:
            if self.plot_hohmann.get():
                if str(planet) == str(self.destination.get(self.destination.curselection())) or \
                   str(planet) == str(self.origin.get(self.origin.curselection())):
                    continue
            planets[planet][4] = np.random.randint(0, 360)  # Random angle between 0 and 360 degree
        for planet, props in planets.items():
            if props[0]:
                x, y = lf.leap_frog(int(nsteps), props[1], props[2], props[4])
                plt.plot(x, y, marker=None, color=props[-1], alpha=0.7)
                try:
                    plt.plot(x[-1], y[-1], marker='o', ms=8, color=props[-1])
                except IndexError:
                    print('Flight Duration cannot be zero.')
                    return 0
            # Sun
            sun_x = 0.
            sun_y = 0.
            shift = -4  # in points
            plt.plot(sun_x, sun_y, marker='o', ms=10, color='Yellow')

        # Plot Hohmann
        if self.plot_hohmann.get():
            print('Fueling Rocket..')
            print('Launching Rocket..')
            transfer_x, transfer_y = hm.hohmann(int(nsteps), planets[self.origin.get(self.origin.curselection())][1], planets[self.origin.get(self.origin.curselection())][2],
                                                planets[self.destination.get(self.destination.curselection())][1], planets[self.destination.get(self.destination.curselection())][2], self.second.get())
            plt.plot(transfer_x, transfer_y, '--', color='Maroon')
        # Set plot limits to largest radius + x
        lim = 0.4
        for values in planets.values():
            if values[0] and values[1] > lim:
                lim = values[1]
        plt.xlabel('x / AU')
        plt.ylabel('y / AU')
        plt.axis("equal")
        plt.xlim(-lim*1.5, lim*1.5)
        plt.ylim(-lim*1.5, lim*1.5)
        plt.grid()
        plt.show()

    def movie(self):
        # Have to select origin and destination of Hohmann Transfer. If not,
        # this function just returns a message and does nothing
        if self.plot_hohmann.get():
            try:
                a, b = self.origin.get(self.origin.curselection()), self.destination.curselection()
            except TclError:
                print('You have to select origin and destiantion for the Hohmann Transfer Orbit')
                return 0

        print('3.. 2.. 1.. Liftoff! ..')

        # Planet: Plotting Boolean, semi-major axisin AU, eccentricity, Offset Angle, color, Position array
        planets = {'Mercury': [self.mercury.get(), 0.387, 0.205, 0.241, 0., 'Gold', []], 'Venus': [self.venus.get(), 0.723, 0.007, 0.615, 0., 'Coral', []],
                   'Earth': [self.earth.get(), 1., 0.017, 1., 0., 'DarkBlue', []], 'Mars': [self.mars.get(), 1.524, 0.094, 1.88, 0., 'Crimson', []],
                   'Jupiter': [self.jupiter.get(), 5.203, 0.049, 11.9, 0., 'orange', []], 'Saturn': [self.saturn.get(), 9.58, 0.057, 29.5, 0., 'Khaki', []],
                   'Uranus': [self.uranus.get(), 19.20, 0.046, 84, 0., 'Turquoise', []], 'Neptune': [self.neptune.get(), 30.06, 0.011, 164.79, 0., 'RoyalBlue', []]}

        # Make save directory if necessary
        save_path = self.save_p.get()
        os.makedirs(save_path, exist_ok=True)

        # Transfer array
        hohmann = []

        # counter for images
        i = 0
        nsteps = 7 * 15 * int(self.duration.get())
        # Set planet angle equal to calculated required angular offset
        if self.plot_hohmann.get():
            planets[self.destination.get(self.destination.curselection())][4] = self.offs.get()
        for planet in planets:
            if self.plot_hohmann.get():
                if str(planet) == str(self.destination.get(self.destination.curselection())) or \
                   str(planet) == str(self.origin.get(self.origin.curselection())):
                    continue
            planets[planet][4] = np.random.randint(0, 360)  # Random angle between 0 and 360 degree

        # If planet is plotted:
        # Create generators to calculate positions of planets
        for planet, properties in planets.items():
            if properties[0]:
                properties.append(all_in_one.leap_frog(nsteps, properties[1], properties[2], properties[4]))

        for x, y in all_in_one.hohmann(nsteps, planets[self.origin.get(self.origin.curselection())][1], planets[self.origin.get(self.origin.curselection())][2],
                                       planets[self.destination.get(self.destination.curselection())][1], planets[self.destination.get(self.destination.curselection())][2], self.second.get()):
            hohmann.append((x, y))
            # If planet is plotted:
            # Calculate next coordinates and append to position array
            for planet, properties in planets.items():
                if properties[0]:
                    properties[6].append(next(properties[7]))
            i += 1
            nth = int(nsteps/100)
            if i % nth == 0 or i == nsteps-1:
                fig = plt.Figure()
                for planet, properties in planets.items():
                    if planets[planet][0]:
                        plt.plot(*zip(*properties[6]), marker=None, color=properties[5], alpha=0.7)
                        plt.plot(*zip(properties[6][-1]), 'o', ms=8, color=properties[5])
                # Sun
                plt.plot(0, 0, 'yo', ms=10)
                # Transfer Orbit
                plt.plot(*zip(*hohmann), '--', color='Maroon')
                plt.plot(*zip(hohmann[-1]), ms=5, color='Maroon')

                # Set plot limits to largest radius + x
                lim = 0.4
                mod = 1.4
                for values in planets.values():
                    if values[0] and values[1] > lim:
                        lim = values[1] * mod
                plt.grid()
                plt.xlabel("x / AU")
                plt.ylabel("y / AU")
                plt.axis("equal")
                plt.axis([-lim, lim, -lim, lim])
                plt.savefig(save_path + str(i) + '.png')
                fig.clear()
                plt.cla()
                plt.clf()

        moviebar = "Shooting movie.. {bar}"
        images = []
        image_files = [save_path + f for f in os.listdir(save_path) if os.path.isfile(os.path.join(save_path, f)) and f[-4:] == '.png']
        image_files.sort()
        for filename in minibar.bar(image_files, template=moviebar):
            images.append(imageio.imread(filename))
            imageio.mimsave(save_path + 'HohmannTransfer.gif', images)
        print('\nDone!')

root = tk.Tk()
root.wm_title('SolarSystemOrbiter')
gui = App(root)
root.mainloop()
