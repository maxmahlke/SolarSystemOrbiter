#! /usr/bin/env python
try:
    import tkinter as tk
    from tkinter import ttk
    from _tkinter import TclError
except ImportError:  # Python 2
    import Tkinter as tk
    from _tkinter import TclError
    import ttk

import calculations
import imageio
import os
import numpy as np
import matplotlib
import collections
matplotlib.use("TkAgg")
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

        self.origin = tk.Listbox(self.transfer_frame, exportselection=0)
        self.origin.grid(row=5, column=0, sticky='EW')
        planets = ['Neptune', 'Uranus', 'Saturn', 'Jupiter', 'Mars', 'Earth', 'Venus', 'Mercury']
        for planet in planets:
            self.origin.insert(0, planet)

        self.destination = tk.Listbox(self.transfer_frame, exportselection=0)
        self.destination.grid(row=5, column=2, sticky='EW')
        for planet in planets:
            self.destination.insert(0, planet)

        self.plot_hohmann = tk.BooleanVar()

        self.plot_hohmann_check = ttk.Checkbutton(self.transfer_frame, text='Plot Hohmann Transfer Orbit', variable=self.plot_hohmann)
        self.plot_hohmann_check.grid(row=6, column=0, sticky='EW', padx=5, pady=5)


        # Mode selection: Orbit Injection, GAM, Nothing
        self.mode = tk.IntVar()

        tk.Radiobutton(self.transfer_frame, text='No Impulse', variable=self.mode, value=0).grid(row=6, column=2, sticky='EW')
        tk.Radiobutton(self.transfer_frame, text='Orbit Injection', variable=self.mode, value=1).grid(row=7, column=0, sticky='EW')
        tk.Radiobutton(self.transfer_frame, text='GAM', variable=self.mode, value=2).grid(row=7, column=2, sticky='EW')

        # Planet offset in degree
        self.offs = tk.DoubleVar()
        # Boolean for acceleration at apohelion
        self.second = tk.BooleanVar()
        # File path for output images
        self.save_p = tk.StringVar()
        self.save_p.set(os.getcwd() + '/MyBigTrip/')
        tk.Entry(self.transfer_frame, textvariable=self.offs).grid(row=8, column=0, sticky='EW')
        tk.Label(self.transfer_frame, text='Planet Offset in Degree').grid(row=8, column=2, sticky='E')
        self.duration = tk.DoubleVar()
        tk.Entry(self.transfer_frame, textvariable=self.duration).grid(row=9, column=0, sticky='EW')
        tk.Label(self.transfer_frame, text='Flight Duration in Earth Weeks').grid(row=9, column=2, sticky='E')

        self.calc = ttk.Button(master, text='Get Rocket Trajectory', command=self.calculate)
        self.calc.grid(column=0, row=10, sticky='EW')

        self.plot = ttk.Button(master, text='Plot', command=self.plot)
        self.plot.grid(column=2, row=10, sticky='EW')

        tk.Entry(self.master, textvariable=self.save_p).grid(row=11, column=0)
        self.movie = ttk.Button(master, text='Make Movie', command=self.movie)
        self.movie.grid(column=2, row=11, sticky='EW')

        # Set-up intital scenario for easy use
        self.mode.set(2)
        self.plot_hohmann.set(True)
        self.earth.set(True)
        self.mars.set(True)
        self.jupiter.set(True)
        self.saturn.set(True)
        self.uranus.set(True)
        self.origin.select_set(2)
        self.destination.select_set(5)
        self.duration.set(800.)
        self.offs.set(105.76)

    def solar_system_planets(self):
       # Planet: Plotting Boolean, semi-major axis in AU, eccentricity, orbital period, angle, color for plot, mass of planet, list of positions (calculated later)
        planets = {'Mercury': {'plotting': self.mercury.get(), 'semi_major': 0.387, 'eccentricity': 0.205, 'orbital_period': 0.241,'angular_offset':  0., 'color': 'Gold', 'mass': 3.302e23, 'position':[]},
                   'Venus': {'plotting': self.venus.get(), 'semi_major': 0.723, 'eccentricity': 0.007, 'orbital_period': 0.615,'angular_offset':  0., 'color': 'Coral', 'mass': 4.8685e24, 'position':[]},
                   'Earth': {'plotting': self.earth.get(), 'semi_major': 1., 'eccentricity': 0.017, 'orbital_period': 1.,'angular_offset':  0., 'color': 'DarkBlue', 'mass': 5.9736e24, 'position':[]},
                   'Mars': {'plotting': self.mars.get(), 'semi_major': 1.524, 'eccentricity': 0.094, 'orbital_period': 1.88,'angular_offset':  0., 'color': 'Crimson', 'mass': 6.4185e23, 'position':[]},
                   'Jupiter': {'plotting': self.jupiter.get(), 'semi_major': 5.203, 'eccentricity': 0.049, 'orbital_period': 11.9,'angular_offset':  0., 'color': 'orange', 'mass': 1.899e27, 'position':[]},
                   'Saturn': {'plotting': self.saturn.get(), 'semi_major': 9.58, 'eccentricity': 0.057, 'orbital_period': 29.5,'angular_offset':  0., 'color': 'Khaki', 'mass': 5.6846e26, 'position':[]},
                   'Uranus': {'plotting': self.uranus.get(), 'semi_major': 19.20, 'eccentricity': 0.046, 'orbital_period': 84, 'angular_offset': 0., 'color': 'Turquoise', 'mass': 8.6832e25, 'position':[]},
                   'Neptune': {'plotting': self.neptune.get(), 'semi_major': 30.06, 'eccentricity': 0.011, 'orbital_period': 164.79,'angular_offset':  0., 'color': 'RoyalBlue', 'mass': 1.0243e26, 'position':[]}}
        return planets # dictionary of dictionaries with planet stats (used namedtuples before, but the attributes have to be mutable)

    def calculate(self):
        # Have to select origin and destination of Hohmann Transfer. If not,
        # this function just returns a message and does nothing
        try:
            a, b = self.origin.get(self.origin.curselection()), self.destination.curselection()
        except TclError:
            print('You have to select origin and destiantion for the Hohmann Transfer Orbit')
            return 0

        print('\nCalculating Rocket trajectory..')

        # Read in planet selection for Hohmann Transfer from UI
        origin_planet = self.origin.get(self.origin.curselection())
        target_planet = self.destination.get(self.destination.curselection())

        # Create Planet quasi-class with function
        planets = self.solar_system_planets()

        # Calculates travel time of rocket and required angular offset of planets
        d_o = planets[origin_planet]['semi_major']
        d_d = planets[target_planet]['semi_major']
        T_d = planets[target_planet]['orbital_period']

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

        # Set values in UI fields for Angular Offset and Flight Duration
        self.duration.set(round(t, 2))
        self.offs.set(round(off, 2))

    def plot(self):
        # Have to select origin and destination of Hohmann Transfer if hohmann should be plotted. If not,
        # this function just returns a message and does nothing
        perform_hohmann_transfer = self.plot_hohmann.get()
        if perform_hohmann_transfer:
            try:
                a, b = self.origin.get(self.origin.curselection()), self.destination.curselection()
            except TclError:
                print('You have to select origin and destiantion for the Hohmann Transfer Orbit')
                return 0

        # Read in planet selection for Hohmann Transfer from UI
        origin_planet = self.origin.get(self.origin.curselection())
        target_planet = self.destination.get(self.destination.curselection())

        # Create Planet quasi-class with function
        planets = self.solar_system_planets()

        # Integration steps. Number of steps is trade-off between accuracy and computational expense
        steps = int(7 * 600 * self.duration.get())

        # Set offset angle of destination planet for HT timing
        # and eccentricity of origin and destinaiton to zero (HT assumes circular orbits)
        if perform_hohmann_transfer:
            planets[target_planet]['angular_offset'] = self.offs.get()
            planets[target_planet]['eccentricity'] = 0.
            planets[origin_planet]['eccentricity'] = 0.

        # Set planet offsets to random values, unless desitnation or target planet
        # Hohmann orbit optional
        for planet in planets:
            if perform_hohmann_transfer:
                if planet in [target_planet, origin_planet]:  # for HT planets, offset already set above
                    continue
            planets[planet]['angular_offset'] = np.random.randint(0, 360)  # Random angle between 0 and 360 degree
        
        # If planet is plotted:
        # Create generators to calculate positions of planets
        # and plot limits
        # Set plot limits to largest radius + x
        lim = 0.4
        for planet in planets:
            if planets[planet]['plotting'] or  planet in [origin_planet, target_planet]:
                if planets[planet]['semi_major'] > lim:
                        lim = planets[planet]['semi_major']
                planets[planet]['generator'] = calculations.planets(steps, planets[planet]['semi_major'], planets[planet]['eccentricity'], planets[planet]['angular_offset'])
        lim *= 1.5

        # Calculate positions of planets and space-craft. Using iterators instead of functions allows to check for current positions and perform GAMs
        
        hohmann = [] # List of spacecraft positions for HT

        # The mode variable defines if we want to do Second Impulse HT, GAM, nothing when arriving at the target planet

        # The routine immediately starts comparing the positions of space-craft and target planet to
        # find out if we have reached the sphere of incfluence. Since we have to update this value, we 
        # use the generator send() function to change it at each integration step
        movie = False
        hohmann_transfer = calculations.hohmann(steps, planets[origin_planet], planets[target_planet], self.mode.get(), movie)
        # initiate generator
        hohmann_transfer.send(None)
        hohmann_transfer.send((0, 0))  # fake target starting position
        # appenrently cannot update generator we iterate over..
        for step in range(0, steps):
            #
            #hohmann.append(next(hohmann_transfer))
        
            # If planet is plotted:
            # Calculate next coordinates and append to position array
            for planet in planets:
                if planets[planet]['plotting'] or planet in [origin_planet, target_planet]:
                    planets[planet]['position'].append(next(planets[planet]['generator']))
            hohmann.append(hohmann_transfer.send(planets[target_planet]['position'][-1]))
            hohmann.append(hohmann_transfer.send(planets[target_planet]['position'][-1]))
        # Remove fake position from target planet
        #planets[target_planet]['position'] = planets[target_planet]['position'][1:]
        
        #clean weird None yields
        hohmann = [h for h in hohmann if h != None]
        # Plot the Solar System simulation
        fig = plt.figure(2)
        for planet, properties in planets.items():
            if planets[planet]['plotting']:
                try:
                    plt.plot(*zip(*planets[planet]['position']), marker=None, color=planets[planet]['color'], alpha=0.7)
                    plt.plot(*zip(planets[planet]['position'][-1]), marker='o', ms=8, color=planets[planet]['color'])
                except IndexError:
                    print('Encounted IndexError when plotting. Have you entered a non-zero, positive flight-duration?')
                    return 0
        # Sun
        plt.plot(0, 0, 'yo', ms=10)
        # Transfer Orbit
        plt.plot(*zip(*hohmann), linestyle='--', color='Maroon')
        plt.plot(*zip(hohmann[-1]), ms=5, color='Maroon')
        plt.grid()
        plt.title('Solar System Simulation')
        plt.xlabel('x / AU')
        plt.ylabel('y / AU')
        plt.axis('equal')
        plt.axis([-lim, lim, -lim, lim])
        print('\nOrbit Transfer\t\t|\t Completed')
        print('Origin Planet\t\t|\t {:s}'.format(origin_planet))
        print('Target Planet\t\t|\t {:s}'.format(target_planet))
        print('Transfer Time\t\t|\t {:.1f} Earth Years'.format(steps / 7 / 150 / 52))
        print('Closest Planet\t\t|\t To Be Implemented')
        # Calculate travelled distance
        travelled_distance = 0
        for i in range(len(hohmann)):
        	try:
        		travelled_distance += np.sqrt((hohmann[i][0]-hohmann[i+1][0])**2 + (hohmann[i][1]-hohmann[i+1][1])**2)
        	except IndexError:
        		break  # reached end of array
        print('Transfer Distance\t|\t {:.1f} AU'.format(travelled_distance))
        fig.show()


    def movie(self):
        # Have to select origin and destination of Hohmann Transfer if hohmann should be plotted. If not,
        # this function just returns a message and does nothing
        perform_hohmann_transfer = self.plot_hohmann.get()
        if perform_hohmann_transfer:
            try:
                a, b = self.origin.get(self.origin.curselection()), self.destination.curselection()
            except TclError:
                print('You have to select origin and destiantion for the Hohmann Transfer Orbit')
                return 0

        # Read in planet selection for Hohmann Transfer from UI
        origin_planet = self.origin.get(self.origin.curselection())
        target_planet = self.destination.get(self.destination.curselection())

        # Create Planet quasi-class with function
        planets = self.solar_system_planets()

        print('3.. 2.. 1.. Liftoff! ..')
        # Make save directory if necessary
        save_path = self.save_p.get()
        try:
            os.makedirs(save_path, exist_ok=True)
        except TypeError:  # Python 2
            if not os.path.exists(save_path):
                os.mkdir(save_path)
        # Transfer array
        hohmann = []

        # counter for images
        i = 0
        # Integration steps. Number of steps is trade-off between accuracy and computational expense        
        steps = 7 * 600 * int(self.duration.get())
        # Set planet angle equal to calculated required angular offset
        if perform_hohmann_transfer:
            planets[target_planet]['angular_offset'] = self.offs.get()
            planets[target_planet]['eccentricity'] = 0.
            planets[origin_planet]['eccentricity'] = 0.

        for planet in planets:
            if perform_hohmann_transfer:
                if planet in [target_planet, origin_planet]:  # for HT planets, offset already set above
                    continue
            planets[planet]['angular_offset'] = np.random.randint(0, 360)  # Random angle between 0 and 360 degree

        # If planet is plotted:
        # Create generators to calculate positions of planets
        # and plot limits
        # Set plot limits to largest radius + x
        lim = 0.4
        for planet in planets:
            if planets[planet]['plotting'] or planet in [target_planet, origin_planet]:
                if planets[planet]['semi_major'] > lim:
                        lim = planets[planet]['semi_major']
                planets[planet]['generator'] = calculations.planets(steps, planets[planet]['semi_major'], planets[planet]['eccentricity'], planets[planet]['angular_offset'])
        lim *= 1.5

        # The routine immediately starts comparing the positions of space-craft and target planet to
        # find out if we have reached the sphere of incfluence. 
        movie = True
        hohmann_transfer = calculations.hohmann(steps, planets[origin_planet], planets[target_planet], self.mode.get(), movie)
        # initiate generator
        hohmann_transfer.send(None)
        hohmann_transfer.send((0, 0))  # fake target starting position
        nth = int(steps/100)
        for step in range(0, steps):

            # If planet is plotted:
            #     Calculate next coordinates and append to position array
            for planet in planets:
                if planets[planet]['plotting'] or planet in [origin_planet, target_planet]:
                    planets[planet]['position'].append(next(planets[planet]['generator']))
            hohmann.append(hohmann_transfer.send(planets[target_planet]['position'][-1]))
            hohmann.append(hohmann_transfer.send(planets[target_planet]['position'][-1]))

            i += 1
            if i % nth == 0 or i == steps-1:
                fig = plt.Figure()
                #clean weird None yields
                hohmann = [h for h in hohmann if h != None]

                for planet, properties in planets.items():
                    if planets[planet]['plotting']:
                        plt.plot(*zip(*planets[planet]['position']), marker=None, color=planets[planet]['color'], alpha=0.7)
                        plt.plot(*zip(planets[planet]['position'][-1]), marker='o', ms=8, color=planets[planet]['color'])
                # Sun
                plt.plot(0, 0, 'yo', ms=10)
                # Transfer Orbit
                plt.plot(*zip(*hohmann), linestyle='--', color='Maroon')
                plt.plot(*zip(hohmann[-1]), ms=5, color='Maroon')
                plt.grid()
                plt.xlabel('x / AU')
                plt.ylabel('y / AU')
                plt.axis('equal')
                plt.axis([-lim, lim, -lim, lim])
                plt.savefig(save_path + str(i) + '.png')
                fig.clear()
                plt.cla()
                plt.clf()
        #clean weird None yields
        hohmann = [h for h in hohmann if h != None]

        images = []
        image_files = [save_path + f for f in os.listdir(save_path) if os.path.isfile(os.path.join(save_path, f)) and f[-4:] == '.png']
        image_files.sort()
        # Skip creating gif for now. Takes to long and messes up image order
        #for filename in image_files:
        #    images.append(imageio.imread(filename))
        #    imageio.mimsave(save_path + 'HohmannTransfer.gif', images)
        print('\nDone!')

root = tk.Tk()
root.wm_title('SolarSystemOrbiter')
gui = App(root)
root.mainloop()
