import time
import math
import numpy as np
import matplotlib.pyplot as plt

class NumericalScheme():
    def __init__(self, title = "", dt=1/100, a=1, N=100, draw=False):
        self.u = np.zeros(N)
        self.dx = 1/N
        self.N = N
        self.dt = dt
        self.a = a
        self.draw = draw
        
        self.i = 0
        
        
        if self.draw:
            fig, ax = plt.subplots(figsize=(8,8))
            
            self.fig = fig
            self.ax = ax
        
            ax.set_title(title)
            ax.set_xlim(0, 1)
            
    def get_t(self):
        return self.i * self.dt
            
    def pass_numerical_scheme(self, scheme, left_k=1, right_k=1, periodic=True):
        self.numerical_scheme = scheme
        self.left_k = left_k
        self.right_k = right_k
        self.periodic = periodic # assume true for now!
    
    def set_initial_condition(self, init_condition):
        for i in range(self.N):
            self.u[i] = init_condition(i / self.N)
            
        self.init_condition = init_condition
        
    def draw(self):
        assert(self.draw, "Need to have draw enabled")
        
        """Draw to the internal figure (for animation)"""
        self.draw_to(self.ax)
        
        plt.ion()
        plt.show()

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
    
    def update(self):
        self.i += 1
        new_u = np.zeros(self.N)
        
        for i in range(self.N):
            data_arr = []
            for j in range(-self.left_k, +self.right_k + 1):
                i_offset = (i + j) % self.N
                
                data_arr.append(self.u[i_offset])
            
            new_u[i] = self.numerical_scheme(data_arr, self.dt, self.dx)
            
        self.u = new_u
            
    
    def draw_to(self, ax):
        """Draw the current state to a given axes object"""
        x = np.arange(0, 1, self.dx)
        
        # Clear the axes first
        ax.clear()
        
        # Draw initial condition
        ax.plot(x, [self.init_condition(xi) if hasattr(self, 'init_condition') else 0 
                    for xi in x], '--', color='gray', alpha=0.5, label='Initial')
        
        # Draw current state
        ax.plot(x, self.u, 'x-', color='purple', label='Current')
        
        ax.set_xlim(0, 1)
        ax.set_xlabel("x")
        ax.set_ylabel("u")
        #ax.set_title(self.ax.get_title())  # Use the title from init
        ax.legend()
        


# numerical schemes

alpha = 0.5

a = 1

def diffuse_data(data, dt, dx):
    u_minus, u, u_plus = data[0], data[1], data[2]
    return u * (1 - alpha) + (u_minus + u_plus) * alpha / 2

def forward_difference(data, dt, dx):
    u, u_plus = data[0], data[1]
    
    return u - a *  (dt / dx) * (u_plus - u)

def backward_difference(data, dt, dx):
    u_minus, u = data[0], data[1]
    
    return u - a * (dt / dx) * (u - u_minus)

def central_difference(data, dt, dx):
    u_minus, u, u_plus = data[0], data[1], data[2]
    
    return u - a * (dt / (2 * dx)) * (u_plus - u_minus)

# initial conditions

def indicator_map(x):
    if 0.25 <= x and x < 0.75:
        return 1

    return 0

def sine_map(x):
    return math.sin(2 * math.pi * x)

#FD_Scheme = NumericalScheme("Backward difference", dt = 1/200)
#FD_Scheme.pass_numerical_scheme(backward_difference)
#FD_Scheme.set_initial_condition(sine_map)


#while True:
#    FD_Scheme.draw()
#    time.sleep(0)
#    FD_Scheme.update()