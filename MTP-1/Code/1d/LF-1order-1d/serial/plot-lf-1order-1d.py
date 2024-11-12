import matplotlib.pyplot as plt

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
        Functiont to read a text file and return it in a list
    """
    data = []

    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            data.append(list(map(float, line.split(' ')[:-1]))[1:-1])

    return data

def plot(density: list[list], grid: list, **kwargs) -> None:
    """
        Funntion plots the initial and final density
    """
    plt.plot(grid, density[0], color="red", label="initial density")
    plt.plot(grid, density[1], color="blue", label="final density")

    plt.title(kwargs["PROBLEM"]+" : "+kwargs["SCHEME"])

    plt.legend()

    plt.savefig(f"{kwargs['PROBLEM']} (DENSITY PLOT).png")

    plt.show()

grid = read("grid.txt")[0]
density = read("density.txt")

params = read_parameters("parameters.txt")

plot(density, grid, **params)

plt.show()



