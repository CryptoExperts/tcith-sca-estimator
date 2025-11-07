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

## MPCitH parameters
ell = 1
slack = 0
masking_order = 0
has_pr_shares = True
has_pr_shares_1seed = False
has_aes_prg = False
accel=True

target_N = 256
N_min = 16

fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()

masking_order_range = [2**i-1 for i in range(0, 6)]

perf = {
    "hash_u":[],
    "hash_m":[],
    "prg_u":[],
    "xof_u":[],
    "sh_eval":[],
    "mult_u":[],
    "mult_m":[],
    #"m_tree":[],
    "refresh":[]
}

timings = []

kappa = 128
ref = 0

print("has_pr_shares", has_pr_shares)
print("has_pr_seed_1seed", has_pr_shares_1seed)
print("accel", accel)

col = 0
for masking_order in masking_order_range:
    slack = (masking_order-1)//2
    slack = 0
    
    x = []
    y = []
    z = []
    best_time = (pow(2, 128), pow(2, 128))
    best_size = (pow(2, 128), pow(2, 128))
    best_size_param = 0
    best_time_param = 0
    #for N in [32]: #256, 512, 1024, 2048]: 
    for N in [N_min + 16*i for i in range(0, 16)] :
    #for N in [target_N]:
        for ell in range(slack+1, slack+2):
            mpc = TCitH_PC_prot(field_size, n, ell)
            mt = TCitH_MT(N, ell, slack, masking_order, tau=None, mpc_protocol=mpc, kappa=kappa, \
                          has_pr_shares=has_pr_shares, has_pr_shares_1seed=has_pr_shares_1seed, has_aes_prg=has_aes_prg, accel=accel)
            tau = mt.parallel_repetitions()
            size = mt.get_size()
            time, detail, keccak_detail = mt.get_performance_index(keccak_perf, aes_perf, isw_perf, refresh_perf, accel=accel)
            p0_u, p0_m, p1, p2_u, p2_m, p3, p4 = detail
            hash_u, hash_m, prg_u, xof_u = keccak_detail

            if masking_order == 0:
                ref = time
            #print(masking_order, time/freq)

            if (N == target_N and ell == slack+1):
                perf["hash_u"].append(100*hash_u/time)
                perf["hash_m"].append(100*hash_m/time)
                perf["prg_u"].append(100*prg_u/time)
                perf["xof_u"].append(100*xof_u/time)
                perf["sh_eval"].append(100*p1/time)
                perf["mult_u"].append(100*p2_u/time)
                perf["mult_m"].append(100*p2_m/time)
                #perf["m_tree"].append(100*p3/time)
                perf["refresh"].append(100*p4/time)  

            if size < best_size[0] :
                best_size = (size, time, time/freq)
                best_size_param = (masking_order, ell, slack, N)
                #print(best_size, best_size_param, tau)
            if time < best_time[1] : 
                best_time = (size, time, time/freq)
                best_time_param = (masking_order, ell, slack, N)
                #print(best_time, best_time_param, tau)
        if N==target_N:
            target_size = size
        if N == target_N:
            timings.append(time)
        x += [best_size[0], best_time[0]]
        y += [best_size[1], best_time[1]]
        z += [best_size[2], best_time[2]]
    ax1.scatter(x,z,label=(masking_order+1), s=dot_size, c=colors_scatter[col], marker=filled_markers[col])
    col += 1
plt.xlabel("Signature size (in bytes)")
ax1.set_ylabel("Time (in seconds)")
#ax1.set_yscale("log")

plt.title("TCitH-MT (MQ-based) masked signature, $N$ from "+str(N_min)+" to "+str(N))
plt.grid(visible=True, which='major', axis='both', linestyle='--', linewidth=.4)
if has_pr_shares == False:
    plt.axvspan(target_size-50, target_size+50, color='grey', alpha=.2)
plt.legend(loc="upper right")
#ax.legend(bbox_to_anchor=(1.1, 1.05))
plt.yscale("log")
plt.savefig("img/tcith-mt-pc-tweak-accel-"+str(N_min)+"-"+str(N)+".pdf", format="pdf", bbox_inches='tight')
plt.show()

#####################################################################

fig, ax = plt.subplots()
width = .6
#fig, ax = plt.subplots()
bottom = np.zeros(len(masking_order_range))
mshares_str = [str(i+1) for i in masking_order_range]

col = 0
for label, num in perf.items():
    if num == [0.0 for i in range(0, len(num))]:
        continue
    if(list(perf.values()).index(num)) < 3:
        p = ax.bar(mshares_str, num, width, label=label, bottom=bottom, color=colors[label], hatch=hatch_mod2[col%2], lw=0, edgecolor='white')
        p = ax.bar(mshares_str, num, width, bottom=bottom, color='none', lw=0, edgecolor='none')
    else:
        p = ax.bar(mshares_str, num, width, label=label, bottom=bottom, color=colors[label], hatch=hatch_mod2[col%2], lw=0, edgecolor='black')
        p = ax.bar(mshares_str, num, width, bottom=bottom, color='none', lw=0, edgecolor='none')
    bottom += num
    #if num != [0.0 for i in range(0, len(num))]:   
    col += 1

plt.xlabel("Masking shares")
plt.ylabel("Performance (in \%)")
plt.title("TCitH-MT (MQ-based) masked signature analysis ($N = "+str(target_N)+"$)")
ax.legend(loc="best")
plt.savefig("img/tcith-mt-pc-tweak-accel-"+str(target_N)+"-analysis.pdf", format="pdf", bbox_inches='tight')
plt.show()