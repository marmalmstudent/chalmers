"""
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
"""
import sys
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
except ImportError:
    sys.stderr.write("Could not execute program: " + str(sys.exc_info()[1])
                     + ".\nPlease install the module and try again")
    sys.stderr.flush()
    sys.exit()
pi = np.pi
c_0 = 2.99792458e8
mu_0 = pi*4e-7
eps_0 = 1/(mu_0*c_0**2)
eta_0 = 377


class Handin5(object):
    def __init__(self):
        """
        Class instantiation
        """

    def simData(self):
        """
        This method updates the data that is to be plotted onto the figure.
        """

    def updatefig(self, simData):
        """
        Handles setting the new image/graph and info text for the simulation.
        Used by the animation function.

        Parameters
        ----------
        simData : tuple
        """
        if (self.taskNbr == 1):
            pass

    def initializeEPlot(self):
        """
        Initializes the 2D plot of the E-field amplitude vs propagated
        distance.
        """
        ymax = 2
        self.figTxt = ax1.text(0, 1.1*ymax, "", color="#000000")
        self.line, = ax1.plot([], [], linestyle="-", color="black")
        ax1.set_ylim(-ymax, ymax)
        ax1.set_xlabel("$z-position [\mu m]$", size=16)
        ax1.set_ylabel("$E-field$", size=16)
        ax1.set_xlim(0, self.dist[np.size(self.dist, axis=0)-1])


if __name__ == "__main__":
    """
    Called when the script runs.
    """
    fig = plt.figure()  # new figure
    hi5 = None
    if (len(sys.argv) > 1):
        args = sys.argv[1:]
        if (args[0] == "1"):
            ax1 = fig.add_subplot(111)
            hi5 = Handin5()  # the simulation class
            hi5.initTask1()  # method for task 1
        else:
            sys.stdout.write("Usage: python <filename.py> <task_nbr>")
            sys.stdout.flush()
            sys.exit()
    else:
        sys.stdout.write("Usage: python <filename.py> <task_nbr>")
        sys.stdout.flush()
        sys.exit()
    # the function handling the the animation itself
    ani = animation.FuncAnimation(fig, hi5.updatefig, hi5.simData,
                                  interval=20, blit=False, repeat=False)
    plt.show()  # plot figure
