import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

def readVector(file_name: str) -> list:
    """
        Function to read a 1d vector from an external file
    """
    data = []
    with open(file_name, mode='r') as f:
        lines = f.readlines()
        for line in lines:
            data.append( list(map(float, line.split(' ')[:-1])) )

    return data[0]

def readResults(file_name: str) -> list[list]:
    """
        Function to read multiple vectors stroed in a single file line by line
    """
    data = []
    with open(file_name, mode='r') as f:
        lines = f.readlines()
        for line in lines:
            data.append( list(map(float, line.split(' ')[:-1]))[1:-1] )

    return data

# Animating
def animate(grid: list, U: list[list]) -> None:
    """
        Function to animate the solution
    """
    def update(i):
        plt.cla()

        plt.plot(grid, U[0], color="red", label="T=0")
        plt.plot(grid, U[i], color="blue", label="T=t")
        plt.xlabel("x")
        plt.ylabel("u(x,t)")

        plt.legend()

    ani = FuncAnimation(plt.gcf(), update, interval=10, frames=len(U), repeat=False)
    ani.save("sine_wave.gif")
    # plt.show()

grid = readVector("grid.txt")
U1 = readResults("u.txt")

animate(grid, U1)