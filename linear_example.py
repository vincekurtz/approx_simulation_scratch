#!/usr/bin/env python

##
#
# A simple demo of using Approximate Bisimulation to find a control
# for a complex system based on a controller designed for a simple 
# system. 
#
##

import numpy as np
import matplotlib.pyplot as plt

class LinearSystem():
    """
    A simple linear system model
    """
    def __init__(self, A, B, C):
        # x' = Ax + Bu
        #  y = Cx
        self.A = np.asarray(A)
        self.B = np.asarray(B)
        self.C = np.asarray(C)

        self.Ts = 0.1  # simulation timestep

    def output(self, x):
        """
        Get the output at state x
        """
        return self.C@x

    def next_state(self, x, u):
        """
        Get the new state resulting from applying control
        u at state x, assuming a sampling time of self.Ts. 
        """
        x_dot = self.A@x + self.B@u
        new_x = x + self.Ts*x_dot
        return new_x

    def simulate(self, x0, u_tape):
        """
        Return the resulting outputs (y) and states (x)
        resulting from applying the control tape u_tape starting
        at initial state x0. 
        """
        states = [x0]
        outputs = [self.output(x0)]

        x = x0
        for t in range(len(u_tape)):
            u = u_tape[t]
            x = self.next_state(x, u)
            y = self.output(x)

            states.append(x)
            outputs.append(y)

        return (states, outputs)

def comparison_plot(output_simple, output_orig):
    """
    Make a plot that compares the outputs of the orignal and the simpliefied systems
    """
    plt.plot(output_simple, label="Reduced Order System")
    plt.plot(output_orig, label="Original System")

    plt.xlabel("time")
    plt.ylabel("output (y)")

    plt.legend()

    plt.show()

def main():
    # define a simple linear system
    sys_simple = LinearSystem([0],[1],[1])   # need to make sure A, B, C are array-like

    # define a double integrator system
    A = np.asarray([[0, 1],[0,0]])
    B = np.asarray([[0],[1]])
    C = np.asarray([1, 0])
    sys_complex = LinearSystem(A,B,C)

    # initial conditions for each system
    x_simp = np.asarray([0])
    x_comp = np.asarray([0,0])

    # simpulate a run with each system
    output_simp = []
    output_comp = []
    for t in range(500):
        output_simp.append(sys_simple.output(x_simp))
        output_comp.append(sys_complex.output(x_comp))

        u_simp = np.asarray([np.sign(np.sin(0.03*t))]) # input to original system 

        # Now we'll calculate the interface to apply an analogous control to the new system
        R = 1
        Q = 0
        P = np.asarray([[1],[0]])
        K = np.asarray([-1, -1])
        u_comp = R*u_simp + Q*x_simp + K@(x_comp - P@x_simp)

        x_simp = sys_simple.next_state(x_simp, u_simp)
        x_comp = sys_complex.next_state(x_comp, u_comp)

    # compare the resulting outputs
    comparison_plot(output_simp, output_comp)
    

if __name__=="__main__":
    main()
