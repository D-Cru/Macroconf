import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa E402
import numpy as np  # noqa E402

time = np.loadtxt(snakemake.input[0], comments=["#", "@"])[:, 0]
potential = np.loadtxt(snakemake.input[0], comments=["#", "@"])[:, 1]

# Plot the potential energy
plt.plot(time, potential)  # , #color='green')
plt.title("Energy minimization. Potential energy over time")
plt.xlabel("Time [ps]")
plt.ylabel("Potential energy [kJ/mol]")
plt.tight_layout()
plt.savefig(snakemake.output[0])
