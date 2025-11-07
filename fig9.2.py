from tcith_estimator_common import *
from tcith_mt_estimator import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

kappa = 128
accel=False
freq = 250*pow(10,6)
dot_size = 30

#matplotlib.use("pgf")
matplotlib.rcParams.update(
    {
        # Adjust to your LaTex-Engine
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "axes.unicode_minus": False,
    }
)
#"""
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

print(CB_color_cycle)

hatch = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']
hatch_mod2 = ['/', '\\']
hatch_mod3 = ['/', '-', '\\']
filled_markers = ['o', 'v', '^', '<', '>', 's', '8', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X', '.']


palette_plasma_x4 = ["#0d0887","#9c179e","#ed7953","#f0f921"]
palette_plasma_x5 = ["#0d0887","#7e03a8","#cc4778","#f89540","#f0f921"]
palette_viridis_x5 = ["#440154","#3b528b","#21918c","#5ec962","#fde725"]
palette_plasma_x6 = ["#0d0887","#6a00a8","#b12a90","#e16462","#fca636","#f0f921"]
palette_plasma_x6_custom = ["#0d0887","#6a00a8","#b12a90","#e16462","#fca636","#fde725"]

palette_plasma_x7 = ["#0d0887","#5c01a6","#9c179e","#cc4778","#ed7953","#fdb42f","#f0f921"]
palette_inferno_x9 = ["#000004","#210c4a","#57106e","#8a226a","#bc3754","#e45a31","#f98e09","#f9cb35","#fcffa4"]
palette_plasma_x9 = ["#0d0887","#4c02a1","#7e03a8","#aa2395","#cc4778","#e66c5c","#f89540","#fdc527","#f0f921"]

colors = {
    "hash_u":"tab:blue",
    "hash_m":"tab:orange",
    "prg_u":"tab:green",
    "prg_m":"tab:red",
    "sh_eval":"yellow",
    "xof_u":"tab:purple",
    "mult_u":"tab:brown",
    "mult_m":"tab:pink",
    "m_tree":"tab:olive",
    "refresh":"tab:cyan"
}

colors = {
    "hash_u" :palette_plasma_x9[0],
    "hash_m" :palette_plasma_x9[1],
    "prg_u"  :palette_plasma_x9[2],
    "prg_m"  :palette_plasma_x9[3],
    "sh_eval":palette_plasma_x9[4],
    "xof_u"  :palette_plasma_x9[5],
    "mult_u" :palette_plasma_x9[6],
    "mult_m" :palette_plasma_x9[7],
    "refresh":palette_plasma_x9[8],
    #"m_tree":"tab:olive",
}


colors_scatter = palette_plasma_x7#["#440154","#414487","#2a788e","#22a884","#7ad151","#fde725"]


#keccak_perf, aes_perf, isw_perf, refresh_perf = m4_perf()
keccak_perf, aes_perf, isw_perf, refresh_perf = riscv_perf(accel)

#MQOM v2 parameters
field_size = 256
N = 256
n = 43

colors_scatter = palette_plasma_x6_custom
from matplotlib.markers import MarkerStyle


## MPCitH parameters
ell = 1
slack = 0


target_N = 32
N_min = 32

max_size = 50000
 
fig, ax1 = plt.subplots()

masking_order=15
kappa = 128
N_range = [N_min + 16*i for i in range(0, 15)]
#ell_range = [(slack+1) + i for i in range(max_slack)]

def generate_points(has_pr_seeds=False,
                    has_pr_seed_1seed=False,
                    has_aes_prg=False,
                    accel=True, slack=0):
    x = []
    y = []
    z = []
    # format best_{size,time}[0, 1, 2] = [size, cycles, freq]
    best_time = (pow(2, 128), pow(2, 128), pow(2, 128))
    best_size = (pow(2, 128), pow(2, 128), pow(2, 128))
    best_size_param = 0
    best_time_param = 0
    for N in N_range :
        for ell in [(slack+1) + i for i in range(7)]:
            mpc = TCitH_PC_prot(field_size, n, ell)
            try:
                mt = TCitH_MT(N, ell, slack, masking_order, tau=None, mpc_protocol=mpc, kappa=kappa, \
                        has_pr_shares=has_pr_seeds, has_pr_shares_1seed=has_pr_seed_1seed, has_aes_prg=has_aes_prg, accel=accel)
            except ValueError:
                continue
            tau = mt.tau
            size = mt.get_size()
            time, _, _ = mt.get_performance_index(keccak_perf, aes_perf, isw_perf, refresh_perf, accel=accel)

            if size < best_size[0] and size < max_size :
                best_size = (size, time, time/freq)
                best_size_param = (ell, slack, N)
                print(best_size, best_size_param, tau)
            if time < best_time[1] and size < max_size : 
                best_time = (size, time, time/freq)
                best_time_param = (ell, slack, N)
                print(best_time, best_time_param, tau)
        if best_size[0] < pow(2, 128) or best_time[0] < pow(2, 128):
            x += [best_size[0], best_time[0]]
            y += [best_size[1], best_time[1]]
            z += [best_size[2], best_time[2]]
    print()
    return x, y, z

x0, y0, z0 = generate_points(False, False, False, True, 0)
x1, y1, z1 = generate_points(True, False, False, True, 0)
x2, y2, z2 = generate_points(False, False, False, True, masking_order)
x3, y3, z3 = generate_points(True, False, False, True, (masking_order-1)//2)

ax1.scatter(x0,z0,label="No tweak", s=dot_size, c=colors_scatter[2], marker=filled_markers[0])
ax1.scatter(x1,z1,label="PR shares", s=dot_size, c=colors_scatter[3], marker=filled_markers[1])
ax1.scatter(x2,z2,label="Slack (full)", s=dot_size, c=colors_scatter[4], marker=filled_markers[2])
ax1.scatter(x3,z3,label="Slack (half) + PR shares", s=dot_size, c=colors_scatter[5], marker=filled_markers[3])

x = x0+x1+x2+x3
z = z0+z1+z2+z3

#ax1.scatter(x,z,label="PR Seeds (d="+str(masking_shares-1)+")", s=dot_size)

x_pareto = []
y_pareto = []

def pareto(x,y,x_ref,y_ref):
    x_pareto = []
    y_pareto = []
    for i in range(len(x)):
        b = True
        for j in range(len(x)):
            for k in range(len(x_pareto)):
                if x[i] == x_pareto[k] and y[i] == y_pareto[k]:
                    b = False
            if x[j] < x[i] and y[j] < y[i]:
                b = False
        if b:
            if x[i] in x_ref and y[i] in y_ref:
                x_pareto.append(x[i])
                y_pareto.append(y[i])
    return x_pareto, y_pareto

x_pareto, z_pareto = pareto(x, z, x, z)

x0_pareto,z0_pareto = pareto(x0,z0,x_pareto,z_pareto)
x1_pareto,z1_pareto = pareto(x1,z1,x_pareto,z_pareto)
x2_pareto,z2_pareto = pareto(x2,z2,x_pareto,z_pareto)
x3_pareto,z3_pareto = pareto(x3,z3,x_pareto,z_pareto)

#print(x_pareto)
#print([round(elem,4) for elem in z_pareto])
ax1.scatter(x0_pareto,z0_pareto, s=dot_size+5, marker=filled_markers[0], alpha=.7, facecolors='none', edgecolors='black')#color='black')
ax1.scatter(x1_pareto,z1_pareto, s=dot_size+5, marker=filled_markers[1], alpha=.7, facecolors='none', edgecolors='black')
ax1.scatter(x2_pareto,z2_pareto, s=dot_size+5, marker=filled_markers[2], alpha=.7, facecolors='none', edgecolors='black')
ax1.scatter(x3_pareto,z3_pareto, s=dot_size+5, marker=filled_markers[3], alpha=.7, facecolors='none', edgecolors='black')
#marker = MarkerStyle('.', fillstyle = 'none')

plt.xlabel("Signature size (in bytes)")
ax1.set_ylabel("Time (in seconds)")
#ax1.set_xlim(left=5000)
#ax1.set_yscale("log")

plt.title("TCitH-MT masked signature, $N$ from $"+str(N_min)+"$ to $"+str(N)+"$, $d="+str(masking_order)+"$")
plt.grid(visible=True, which='major', axis='both', linestyle='--', linewidth=.4)
#plt.axvspan(target_size-100, target_size+100, color='grey', alpha=.2)
plt.legend(loc="upper right")
#ax.legend(bbox_to_anchor=(1.1, 1.05))
plt.yscale("log")
plt.savefig("img/tcith-mt-pc-comparison-"+str(masking_order+1)+"-shares.pdf", format="pdf", bbox_inches='tight')
plt.show()

#####################################################################