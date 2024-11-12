import matplotlib.pyplot as plt

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

def plot(density: list[list], grid: list) -> None:
    """
        Funntion plots the initial and final density
    """
    plt.plot(grid, density[0], color="red", label="initial density")
    plt.plot(grid, density[1], color="blue", label="final density")

    TITLE = "MCW - LW Scheme 2nd Order Serial - 1D"
    # TITLE = "SCW - LW Scheme 2nd Order Serial - 1D"
    # TITLE = "BLW - LW Scheme 2nd Order Serial - 1D"

    plt.title(TITLE)

    plt.legend()

    plt.savefig(f"{TITLE} (DENSITY PLOT).png")

    plt.show()

if __name__ == "__main__":
    grid = read("grid.txt")[0]
    density = read("density.txt")

    plot(density, grid)


