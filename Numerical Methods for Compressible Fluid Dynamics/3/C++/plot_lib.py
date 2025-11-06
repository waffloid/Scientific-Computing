"""
goals: 
    * plot solution given .dat
    * animate solution given .dat
"""
from enum import Enum

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import numpy as np
import os
import time

plt.style.use('fivethirtyeight')

Y_AXIS_BUFFER = 0.1

def load_dat(filename):
    data = np.loadtxt(filename)
    return np.array(data)

class Content(Enum):
    SOLUTION = 1
    ERROR = 2
    ANALYTIC = 3

class simulation():
    static_axes = True

    def __init__(self, filename: str):
        self.init_condition = load_dat(filename + "/" + "initial_condition.dat") 
        self.N = len(self.init_condition)
        self.sim_length = len(os.listdir(filename + '/simulation'))
        self.filename = filename

        self.simulation = np.zeros([self.sim_length, self.N])
        self.error = np.zeros([self.sim_length, self.N])

        for i in range(self.sim_length):
            self.simulation[i] = load_dat(filename + "/simulation/" + str(i) + ".dat")
            self.error[i] = np.zeros(len(self.simulation[i]))# load_dat(filename + "/error/" + str(i) + ".dat")
                                     

    def draw_to(self, frame: int, ax, content=Content.SOLUTION):
        source = []
        data = []

        match content:
            case Content.SOLUTION:
                source = self.simulation
            case Content.ERROR:
                source = self.error
            case Content.ANALYTIC:
                source = self.simulation + self.error
    
        data = source[frame]

        """Draw a given state to a given axes object"""
        x = np.arange(0, self.N)
        
        # Clear the axes first
        ax.clear()
        
        # Draw initial condition
        ax.plot(x / self.N, [self.init_condition[i] if hasattr(self, 'init_condition') else 0 
                    for i in x], '--', color='gray', alpha=0.5, label='Initial')
        
        # Draw current state
        ax.plot(x / self.N, data, 'x', color='purple', label='Current')
        
        ax.set_xlim(0, 1)
        ax.set_xlabel("x")
        ax.set_ylabel("u")

        if self.static_axes:
            all_min = source.min()
            all_max = source.max()
            data_range = all_max - all_min

            ax.set_ylim(all_min - Y_AXIS_BUFFER * data_range, 
                        all_max + Y_AXIS_BUFFER * data_range)

            print(data_range, all_min, all_max)
        
        ax.legend()

    def animate(self, content = Content.SOLUTION, interval=100):
        fig, ax = plt.subplots()
        
        def update(frame):
            self.draw_to(frame, ax, content)

        anim = FuncAnimation(fig, update, frames=self.sim_length, interval=interval, repeat=True)

        plt.show()
        return anim
        
    def plot_l1_error(self):
        fig, ax = plt.subplots()

        ax.set_title("Error")
        ax.set_xlabel("N")
        ax.set_ylabel("L1 error")

        n = np.arange(0, self.N)
        errors = np.zeros(self.N)

        for i in range(self.N):
            print("Loading " + str(i))
            errors[i] = load_dat(self.filename + '/error_l1/' + str(i) + ".dat")

        plt.plot(n, errors)
        ax.legend()
        plt.show()




FD_advection = simulation("Output/Burger FORCE")


#FD_advection.plot_l1_error()
questionable_anim = FD_advection.animate(content=Content.SOLUTION)



#fig, ax = plt.subplots(figsize=(8,8))
#FD_advection.draw_to(10, ax)
#plt.show()

#todo: get simulation length from the data, make sure that cpp code actually compiles and creates dirs, make sure sim_lnegth & simulation etc are actually properties in the class


