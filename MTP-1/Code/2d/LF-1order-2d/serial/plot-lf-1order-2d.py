import matplotlib.pyplot as plt
import numpy as np

def read_parameters(filename: str) -> dict:
    """
        Function to get the parameters for the current problem
    """
    params = {}
    with open(filename, "r") as f:
        lines = f.readlines()
        params[f"N"] = int(lines[0][:-1])
        params[f"PROBLEM"] = lines[1][:-1]
        params[f"BC"] = lines[2][:-1]
        params[f"TIME"] = list(map(float, lines[3][:-1].split(" ")))
        params[f"SCHEME"] = lines[4][:-1]

    return params

def read(filename: str) -> list:
    """
        Function to read a text file and return it in a list
    """
    data = []

    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            data.append(list(map(float, line.split(' ')[:-1])))

    return data

def plot2d(density0: list[list], density1: list[list], gridX: list, gridY: list, **kwargs) -> None:
    """
        Function to make a 2D plot
    """
    plt.imshow(density0, extent=[gridX[0], gridX[-1], gridY[0], gridY[-1]], label=f"density (t={kwargs['TIME'][0]})")
    # plt.imshow(density1, extent=[gridX[0], gridX[-1], gridY[0], gridY[-1]], label=f"density (t={kwargs['TIME'][1]})")
    plt.title(kwargs["PROBLEM"] + " : " + kwargs["SCHEME"])
    
    plt.show()

def plot3d(density: list[list], gridX: list, gridY: list, **kwargs) -> None:
    """
        Function to make a 3D surface plot
    """
    X, Y = np.meshgrid(gridX, gridY)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, np.array(density))

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    ax.set_title(kwargs["PROBLEM"] + " : " + kwargs["SCHEME"])

    plt.show()

gridX = read("grid.txt")[0]
gridY = read("grid.txt")[1]
density0 = read("initial_density.txt")
density1 = read("final_density.txt")
params = read_parameters("parameters.txt")

plot2d(density0, density1, gridX, gridY, **params)
plot2d(density1, density1, gridX, gridY, **params)

# plot3d(density, gridX, gridY, **params)
