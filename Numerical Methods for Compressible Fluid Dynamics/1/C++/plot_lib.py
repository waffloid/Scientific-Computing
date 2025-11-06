"""
goals: 
    * plot solution given .dat
    * animate solution given .dat
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import numpy as np
import os
import time

Y_AXIS_BUFFER = 0.1

plt.style.use('fivethirtyeight')

def load_dat(filename):
    data = np.loadtxt(filename)
    return np.array(data)

class simulation():
    static_axes = True

    def __init__(self, filename):
        self.init_condition = load_dat(filename + "/" + "initial_condition.dat") 
        self.N = len(self.init_condition)
        self.sim_length = len(os.listdir(filename + '/simulation'))

        self.simulation = np.zeros([self.sim_length, self.N])

        for i in range(self.sim_length):
            self.simulation[i] = load_dat(filename + "/simulation/" + str(i) + ".dat")

    def draw_to(self, simulation_slice, ax):
        
        """Draw a given state to a given axes object"""
        x = np.arange(0, self.N)
        
        # Clear the axes first
        ax.clear()
        
        # Draw initial condition
        ax.plot(x / self.N, [self.init_condition[i] if hasattr(self, 'init_condition') else 0 
                    for i in x], '--', color='gray', alpha=0.5, label='Initial')
        
        # Draw current state
        ax.plot(x / self.N, self.simulation[simulation_slice], 'x', color='purple', label='Current')
        
        ax.set_xlim(0, 1)
        ax.set_xlabel("x")
        ax.set_ylabel("u")

        if self.static_axes:
            all_min = self.simulation.min()
            all_max = self.simulation.max()
            data_range = all_max - all_min

            ax.set_ylim(all_min - Y_AXIS_BUFFER * data_range, 
                        all_max + Y_AXIS_BUFFER * data_range)

        
        ax.legend()

    def animate(self, interval=100):
        fig, ax = plt.subplots()
        
        def update(frame):
            self.draw_to(frame, ax)

        anim = FuncAnimation(fig, update, frames=self.sim_length, interval=interval, repeat=True)

        plt.show()
        return anim
        
    def plot_error(self):
        fig, ax = plt.subplots()

        ax.set_title("Error")
        ax.set_xlabel("N")
        ax.set_ylabel("L1 error")

        n = np.arange(0, self.N)
        errors = np.zeros(self.N)

        for i in range(self.N):
            errors[i] = load_dat(filename + '/metadata/error_L1_' + i)

        plt.plot(n, errors)
        ax.legend()
        ax.show()




FD_advection = simulation("Lax Friedrichs")
FD_advection.static_axes = False

FD_advection.plot_error()

#questionable_anim = FD_advection.animate()



#fig, ax = plt.subplots(figsize=(8,8))
#FD_advection.draw_to(10, ax)
#plt.show()

#todo: get simulation length from the data, make sure that cpp code actually compiles and creates dirs, make sure sim_lnegth & simulation etc are actually properties in the class


