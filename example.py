from tcith_estimator_common import *
from tcith_mt_estimator import *

### Platform Configuration
# Platform hash acceleration
accel=True
# Platform CPU frequency
freq = 250*pow(10,6)

### Security Parameter
kappa = 128

### SCA Protection
# Masking order
masking_order = 7

### Tweak Configuration
# Slack parameter (0 = no slack)
#slack = (masking_order-1)//2 # half slack
#slack = masking_order # full slack
slack = 0
# Pseudo random shares tweak
has_pr_shares = True

### MPC Protocol
field_size = 256
n = 43 # Input size
N = field_size # Number of parties
ell = 1 # Privacy Threshold. It is required that ell > slack. If slack == 0, we then have ell = 1.
# Instantiate the MPC protocol Pi_PC from FR23 article
mpc = TCitH_PC_prot(field_size, n, ell=ell)
mt = TCitH_MT(N, ell, slack, masking_order, tau=None, mpc_protocol=mpc, kappa=kappa, \
              has_pr_shares=has_pr_shares, accel=accel)
# Given the TCitH instance, compute the required parallel repetitions required to reach kappa bits of security. This step IS REQUIRED in order to reach the target security.
#tau = mt.parallel_repetitions()
# Now that the instance is fully initialized, compute the signature size
size = mt.get_size()
# Obtain the benchmarks for the different functions. See the perf.py file for the benchmakrs format
keccak_perf, aes_perf, isw_perf, refresh_perf = riscv_perf(accel)
# Compute the performance in clock cycles, see the tcith_mt_estimator.py file for the format of the detailed performances. The time variable contains the full performance in clock cycles.
time, detail, keccak_detail = mt.get_performance_index(keccak_perf, aes_perf, isw_perf, refresh_perf, accel=accel)

print(f"signature estimated time: {time/freq:.2f} s")
print(f"signature estimate size: {size} B")