#! /usr/bin/env python

import matplotlib.pyplot as plt
import leap_frog as lf
import tkinter as tk
from tkinter import ttk
from _tkinter import TclError
import hohmann as hm
import imageio
import os
import all_in_one
import numpy as np
import minibar

class app:
    def __init__(self, master):
        self.master = master
        self.master.title = 'Hohmann-Transfer Calculator'

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
        self.origin.grid(row=4, column=0, sticky='EW')
        for item in ['Neptune', 'Uranus', 'Saturn', 'Jupiter', 'Mars', 'Earth', 'Venus', 'Mercury']:
            self.origin.insert(0, item)

        self.destination = tk.Listbox(self.transfer_frame, exportselection=0)
        self.destination.grid(row=4, column=2, sticky='EW')
        for item in ['Neptune', 'Uranus', 'Saturn', 'Jupiter', 'Mars', 'Earth', 'Venus', 'Mercury']:
            self.destination.insert(0, item)

        # Planet offset in degree
        self.offs = tk.DoubleVar()
        # Boolean for acceleration at apohelion
        self.second = tk.BooleanVar()
        # File path for output images
        self.save_p = tk.StringVar()
        self.save_p.set(os.getcwd())

        tk.Entry(self.transfer_frame, textvariable=self.offs).grid(row=6, column=0, sticky='EW')
        tk.Label(self.transfer_frame, text='Planet Offset in Degree').grid(row=6, column=2, sticky='E')
        self.duration = tk.DoubleVar()
        tk.Entry(self.transfer_frame, textvariable=self.duration).grid(row=7, column=0, sticky='EW')
        tk.Label(self.transfer_frame, text='Flight Duration in Earth Weeks').grid(row=7, column=2, sticky='E')

        tk.Checkbutton(self.transfer_frame, text='Acceleration at Apohelion', variable=self.second).grid(row=5, column=0, sticky='EW')

        self.calc = ttk.Button(master, text='Get Rocket Trajectory', command=self.calc)
        self.calc.grid(column=0, row=6, sticky='EW')

        self.plot = ttk.Button(master, text='Plot', command=self.plot)
        self.plot.grid(column=2, row=6, sticky='EW')

        tk.Entry(self.master, textvariable=self.save_p).grid(row=7, column=0)
        self.movie = ttk.Button(master, text='Make Movie', command=self.movie)
        self.movie.grid(column=2, row=7, sticky='EW')

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
        # Using Keplers third law units of 1 AU and 1 Earth year: T^2 = a^3
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
        # Have to select origin and destination of Hohmann Transfer. If not,
        # this function just returns a message and does nothing
        try:
            a, b = self.origin.get(self.origin.curselection()), self.destination.curselection()
        except TclError:
            print('You have to select origin and destiantion for the Hohmann Transfer Orbit')
            return 0

        # Planet: Plotting Boolean, distances to sun in AU, angle
        planets = {'Mercury': [self.mercury.get(), 0.387, 0.241, 0.], 'Venus': [self.venus.get(), 0.723, 0.615, 0.],
                   'Earth': [self.earth.get(), 1., 1., 0.], 'Mars': [self.mars.get(), 1.524,  1.88, 0.],
                   'Jupiter': [self.jupiter.get(), 5.203, 11.9, 0.], 'Saturn': [self.saturn.get(), 9.58, 29.5, 0.],
                   'Uranus': [self.uranus.get(), 19.20, 84, 0.], 'Neptune': [self.neptune.get(), 30.06, 164.79, 0.]}

        nsteps = 7 * 15 * self.duration.get()
        planets[self.destination.get(self.destination.curselection())][3] = self.offs.get()

        for planet, props in planets.items():
            if props[0]:
                x, y = lf.leap_frog(int(nsteps), props[1], props[3])
                plt.plot(x, y)
                plt.plot(x[-1], y[-1], 'o', ms=8)
            # Sun
            plt.plot(0, 0, 'yo', ms=10)

        # Plot Hohmann
        print('Fueling Rocket..')
        print('Launching Rocket..')
        transfer_x, transfer_y = hm.hohmann(int(nsteps), planets[self.origin.get(self.origin.curselection())][1],
                                            planets[self.destination.get(self.destination.curselection())][1], self.second.get())
        plt.plot(transfer_x, transfer_y)
        # Set plot limits to largest radius + x
        lim = 0.4
        for values in planets.values():
            if values[0] and values[1] > lim:
                lim = values[1]
        plt.xlim(-lim*1.1, lim*1.1)
        plt.ylim(-lim*1.1, lim*1.1)
        plt.xlabel('x / AU')
        plt.ylabel('y / AU')
        plt.grid()
        plt.show()

    def movie(self):
        # Have to select origin and destination of Hohmann Transfer. If not,
        # this function just returns a message and does nothing
        try:
            a, b = self.origin.get(self.origin.curselection()), self.destination.curselection()
        except TclError:
            print('You have to select origin and destiantion for the Hohmann Transfer Orbit')
            return 0

        print('Here we go! ..')
        planets = {'Mercury': [self.mercury.get(), 0.387, 0.241, 0.], 'Venus': [self.venus.get(), 0.723, 0.615, 0.],
                   'Earth': [self.earth.get(), 1., 1., 0.], 'Mars': [self.mars.get(), 1.524,  1.88, 0.],
                   'Jupiter': [self.jupiter.get(), 5.203, 11.9, 0.], 'Saturn': [self.saturn.get(), 9.58, 29.5, 0.],
                   'Uranus': [self.uranus.get(), 19.20, 84, 0.], 'Neptune': [self.neptune.get(), 30.06, 164.79, 0.]}

        # Make save directory if necessary
        save_path = self.save_p.get()
        os.makedirs(save_path, exist_ok=True)

        # Planet and Transfer arrays
        hohmann_x = []
        hohmann_y = []
        merx = []
        mery = []
        venx = []
        veny = []
        earx = []
        eary = []
        marx = []
        mary = []
        jupx = []
        jupy = []
        satx = []
        saty = []
        urax = []
        uray = []
        nepx = []
        nepy = []

        # counter for images
        i = 0
        nsteps = 7 * 15 * int(self.duration.get())
        # Set planet angle equal to calculated required angular offset
        planets[self.destination.get(self.destination.curselection())][3] = self.offs.get()

        if planets['Mercury'][0]:
            mercury = all_in_one.leap_frog(nsteps, planets['Mercury'][1], planets['Mercury'][3])
        if planets['Venus'][0]:
            venus = all_in_one.leap_frog(nsteps, planets['Venus'][1], planets['Venus'][3])
        if planets['Earth'][0]:
            earth = all_in_one.leap_frog(nsteps, planets['Earth'][1], planets['Earth'][3])
        if planets['Mars'][0]:
            mars = all_in_one.leap_frog(nsteps, planets['Mars'][1], planets['Mars'][3])
        if planets['Jupiter'][0]:
            jupiter = all_in_one.leap_frog(nsteps, planets['Jupiter'][1], planets['Jupiter'][3])
        if planets['Saturn'][0]:
            saturn = all_in_one.leap_frog(nsteps, planets['Saturn'][1], planets['Saturn'][3])
        if planets['Uranus'][0]:
            uranus = all_in_one.leap_frog(nsteps, planets['Uranus'][1], planets['Uranus'][3])
        if planets['Neptune'][0]:
            neptune = all_in_one.leap_frog(nsteps, planets['Neptune'][1], planets['Neptune'][3])

        for x, y in all_in_one.hohmann(nsteps, planets[self.origin.get(self.origin.curselection())][1],
                                       planets[self.destination.get(self.destination.curselection())][1], self.second.get()):
            hohmann_x.append(x)
            hohmann_y.append(y)

            if planets['Mercury'][0]:
                mx, my = next(mercury)
                merx.append(mx)
                mery.append(my)
            if planets['Venus'][0]:
                vx, vy = next(venus)
                venx.append(vx)
                veny.append(vy)
            if planets['Earth'][0]:
                ex, ey = next(earth)
                earx.append(ex)
                eary.append(ey)
            if planets['Mars'][0]:
                max, may = next(mars)
                marx.append(max)
                mary.append(may)
            if planets['Jupiter'][0]:
                jx, jy = next(jupiter)
                jupx.append(jx)
                jupy.append(jy)
            if planets['Saturn'][0]:
                sx, sy = next(saturn)
                satx.append(sx)
                saty.append(sy)
            if planets['Uranus'][0]:
                ux, uy = next(uranus)
                urax.append(ux)
                uray.append(uy)
            if planets['Neptune'][0]:
                nx, ny = next(neptune)
                nepx.append(nx)
                nepy.append(ny)

            i += 1
            nth = int(nsteps/50)

            if i % nth == 0 or i == nsteps-1:
                fig = plt.Figure()

                if planets['Mercury'][0]:
                    plt.plot(merx, mery, 'b-')
                    plt.plot(mx, my, 'o', ms=8)
                if planets['Venus'][0]:
                    plt.plot(venx, veny, 'b-')
                    plt.plot(vx, vy, 'o', ms=8)
                if planets['Earth'][0]:
                    plt.plot(earx, eary, 'b-')
                    plt.plot(ex, ey, 'o', ms=8)
                if planets['Mars'][0]:
                    plt.plot(marx, mary, 'b-')
                    plt.plot(max, may, 'o', ms=8)
                if planets['Jupiter'][0]:
                    plt.plot(jupx, jupy, 'b-')
                    plt.plot(jx, jy, 'o', ms=8)
                if planets['Saturn'][0]:
                    plt.plot(satx, saty, 'b-')
                    plt.plot(sx, sy, 'o', ms=8)
                if planets['Uranus'][0]:
                    plt.plot(urax, uray, 'b-')
                    plt.plot(ux, uy, 'o', ms=8)
                if planets['Neptune'][0]:
                    plt.plot(nepx, nepy, 'b-')
                    plt.plot(nx, ny, 'o', ms=8)

                # Sun
                plt.plot(0, 0, 'yo', ms=10)
                plt.plot(hohmann_x, hohmann_y, 'r-')
                plt.plot(x, y, 'x', ms=5)

                # Set plot limits to largest radius + x
                plt.grid()
                lim = 0.4
                for values in planets.values():
                    if values[0] and values[1] > lim:
                        lim = values[1]
                plt.xlim(-lim*1.1, lim*1.1)
                plt.ylim(-lim*1.1, lim*1.1)
                plt.xlabel("x / AU")
                plt.ylabel("y / AU")
                plt.savefig(save_path + str(i) + '.png')
                fig.clear()
                plt.cla()
                plt.clf()

        moviebar = "Shooting movie.. {bar}  {eta}"
        images = []
        image_files = [save_path + f for f in os.listdir(save_path) if os.path.isfile(os.path.join(save_path, f)) and f[-4:] == '.png']
        for filename in minibar.bar(image_files, template=moviebar):
            images.append(imageio.imread(filename))
            imageio.mimsave(save_path + 'HohmannTransfer.gif', images)  # Save object dictionary entry to text file
        print('Done!')

root = tk.Tk()
gui = app(root)
root.mainloop()
