import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa E402
import numpy as np  # noqa E402

time = np.loadtxt(snakemake.input[0], comments=["#", "@"])[:, 0]
potential = np.loadtxt(snakemake.input[0], comments=["#", "@"])[:, 1]
temperature = np.loadtxt(snakemake.input[0], comments=["#", "@"])[:, 2]
density = np.loadtxt(snakemake.input[0], comments=["#", "@"])[:, 3]

# Plot the potential energy
plt.figure()
plt.plot(time, potential)  # , #color='green')
plt.title("Equilibration. Potential energy over time")
plt.xlabel("Time [ps]")
plt.ylabel("Potential energy [kJ/mol]")
plt.tight_layout()
plt.savefig(snakemake.output.pot)

plt.figure()
plt.plot(time, temperature)  # , #color='green')
plt.title("Equilibration. Temperature over time")
plt.xlabel("Time [ps]")
plt.ylabel("Temperature [K]")
plt.tight_layout()
plt.savefig(snakemake.output.temp)

plt.figure()
plt.plot(time, density)  # , #color='green')
plt.title("Equilibration. Density over time")
plt.xlabel("Time [ps]")
plt.ylabel("Density [kg/m^3]")
plt.tight_layout()
plt.savefig(snakemake.output.dens)
