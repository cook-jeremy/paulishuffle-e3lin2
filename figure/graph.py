import csv, pickle
import matplotlib.pyplot as plt
import numpy as np

gamma = []
#gamma_offset = 0.025
nest = []
nest_err = []
pauli = []
pauli_hoeff = []
diff = []

with open('cvsgamma_data.csv') as cvsgamma_data_file:
    csv_reader = csv.reader(cvsgamma_data_file, delimiter=",")
    b = 0
    for row in csv_reader:
        if b == 0 or b == 1:
            b += 1
            continue
        row = [float(s) for s in row]
        gamma.append(row[0]/np.pi)
        nest.append(row[1])
        nest_err.append(row[2])
        pauli.append(row[3])
        pauli_hoeff.append(row[4])
        diff.append(abs(row[1] - row[3]))

fig, axs = plt.subplots(2,1)

#def error_plot(x, y, offset, error):
#    ax.errorbar([a + offset for a in x], y, yerr=error, linestyle="None", elinewidth=1, capsize = 1.5, ecolor='black')
#error_plot(gamma, nest, gamma_offset/2, nest_err)
#error_plot(gamma, pauli, -gamma_offset/2, pauli_err)

axs[0].errorbar(gamma, diff, yerr=nest_err, linestyle="None", c="k",  marker="x", label="$|\\langle C \\rangle_{Nest} - \\langle C \\rangle_{Heis}| \pm \\varepsilon_{Nest}$ ")
axs[0].plot(gamma, diff, linestyle="--", c="k")
axs[0].scatter(gamma, pauli_hoeff, marker='o', c="k", label="$\\varepsilon_{Heis}$")
axs[0].plot(gamma, pauli_hoeff, c="k")

axs[0].legend()
axs[0].set_xlabel("$\gamma/\pi$")
axs[0].set_ylabel("Error")
axs[0].grid(b=True, axis='both', linestyle='--')
axs[0].set_xticks(np.arange(0,1.1,step=0.1))


m = []
nest = []
e_nest = []
pauli = []
e_pauli = []
diff = []

with open('30_bounds.csv') as bounds_file:
    csv_reader = csv.reader(bounds_file, delimiter=",")
    b = 0
    for row in csv_reader:
        if b == 0: 
            b = 1
            continue
        row = [float(s) for s in row]
        m.append(10*row[0])
        nest.append(row[1])
        e_nest.append(row[2])
        pauli.append(row[3])
        e_pauli.append(row[4])
        diff.append(abs(row[1] - row[3]))
        
axs[1].errorbar(m, diff, yerr=e_nest, linestyle="None", c="k",  marker="x")
axs[1].plot(m, diff, linestyle="--", c="k")
axs[1].scatter(m, e_pauli, marker='o', c="k")
axs[1].plot(m, e_pauli, c="k")

axs[1].set_xlabel("$m$ = number of equations ")
axs[1].set_ylabel("Error")
axs[1].grid(b=True, axis='both', linestyle='--')
axs[1].set_xticks(np.arange(40,81,step=10))
axs[1].set_yscale("log")

plt.suptitle("Accuracy of Heisenberg propagation")
plt.show()
