#!python3
import numpy as np
import matplotlib.pyplot as plt

markers = ['o', 'd', 's', '+', '^', '<', '>']
n_steps = 200
labels = [f"n = {j+1}" for j in range(len(markers))]

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

def make_legend_handles(ax, color):
    handles = []
    for m in markers:
        h, = ax.plot([], [], marker=m, ls='', mfc="none", mec='k')
        handles.append(h)
    return handles

def plot_history(ax, prefix, color):
    # storage for each marker group
    x_list = [[] for _ in markers]
    y_list = [[] for _ in markers]

    # sum curve
    sum_x, sum_y = [], []

    for i in range(n_steps):
        data = np.load(f"data/A_b/{prefix}_{i}_linear_history.npy")

        for j, val in enumerate(data):
            x_list[j].append(i)
            y_list[j].append(val)

        sum_x.append(i)
        sum_y.append(np.sum(data))

    # plot markers (no lines)
    for j, m in enumerate(markers):
        if x_list[j]:
            ax.plot(x_list[j], y_list[j],
                    marker=m, ls='',   # <--- no line
                    mfc="none", mec='k')

    # twin y-axis: summed curve
    ax_top = ax.twinx()
    ax_top.plot(sum_x, sum_y, color='k', lw=2)
    ax_top.set_ylim([0, 60])
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_ylabel("Total Number of Linear Iteration")

# Temperature
axes[0].set_xlim([-1, 40])
axes[0].set_ylim([0, 25])

plot_history(axes[0], 'cgt', 'red')
axes[0].set_title("Temperature Formulation")
axes[0].legend(make_legend_handles(axes[0], 'red'),
               labels, loc="upper left", ncol=4)

# Conservative
axes[1].set_xlim([-1, 40])
axes[1].set_ylim([0, 25])
axes[0].plot([25, 25], [0,25],'--r')
axes[1].plot([25, 25], [0,25],'--r')
axes[0].plot([20, 20], [0,25],'--b')
axes[1].plot([20, 20], [0,25],'--b')
plot_history(axes[1], 'cgc', 'blue')
axes[1].set_title("Conservative Formulation")
axes[1].legend(make_legend_handles(axes[1], 'blue'),
               labels, loc="upper left", ncol=4)
axes[0].set_xlabel("time step")
axes[1].set_xlabel("time step")
axes[0].set_ylabel("Linear Iterations")

plt.tight_layout()
plt.savefig("linear_iterations_per_timestep_per_formulation.png", dpi=300)
plt.show()
