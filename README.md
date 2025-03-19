Masked TCitH signature performance and size estimator
=====================================================

In this repo we present the performance estimator following the methodology presented in the "Masking-Friendly Post-Quantum Signatures in the Threshold-Computation-in-the-Head Framework" article.

The estimator allows to estimate the performance of a masked signature given a set of parameters and benchmarks of some individual primitives. It only requires `python3` to be executed. In order to get an estimation, the program requires:
- A MPC protocol instance,
- Benchmarks of some primitives ((masked) Keccak, (masked) field multiplication, mask refreshing) at a given masking order,
- A MPCitH instance, instantiated using the MPC protocol instance.

### How to use the estimator?

The estimator is divided in different files, namely:
- `tcith_estimator_common.py`, which contains the parent class for MPCitH, and the different classes for MPC protocols,
- `tcith_ggm_estimator.py`, which contains the class for TCitH-GGM, implementing the parallel trees and slack tweaks as described in the reference paper,
- `tcith_mt_estimator.py`, which contains the class for TCitH-MT, implementing the pseudorandom shares and slack tweaks as described in the reference paper,
- `perf.py`, which contains all the functions concerning the performance benchmarks. In particular, it contains the function `riscv_perf` that returns performance in clock cycles for Keccak, AES (which is not supported at the time), PINI multiplication and mask refreshing for `i = 1..32` masking shares.

An estimator instance can be run as:
```python
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

mt = TCitH_MT(N, ell, slack, masking_order, mpc_protocol=mpc, \
                kappa=kappa, has_pr_shares=has_pr_shares, accel=accel)
# Now that the instance is fully initialized, compute the signature size (in bytes)
size = mt.get_size()
# Obtain the benchmarks for the different functions. See the perf.py file for the benchmakrs format
keccak_perf, aes_perf, isw_perf, refresh_perf = riscv_perf(accel)
# Compute the performance in clock cycles, see the tcith_mt_estimator.py file for the format of the detailed performances. The time variable contains the full performance in clock cycles.
time, detail, keccak_detail = mt.get_performance_index(keccak_perf, aes_perf, isw_perf, refresh_perf, accel=accel)

print(f"signature estimated time: {time/freq:.2f} s")
print(f"signature estimate size: {size} B")

```

### Note on hardware acceleration.

The script allows for estimating the performance using a Keccak hardware accelerator, as described in [Saa24]. The benchmarks for this accelerator are made differently compared with the other primitives. Indeed, while the unaccelerated benchmarks only observe the running time of a given function, the accelerated benchmarks take more parameters into account, such as the size of the input and output of the function (since using the accelerator require memory transfers, they NEED to be taken into account as they are much slower than running the accelerator itself). As a consequence, the benchmarks are made differently. Since the hardware accelerator is used to benchmark SHA3 hash function and XOF, the input and output sizes can not be predicted thus we would need to make many benchmarks to take all the different possibilities into account, and it would be sufficient yet for an arbitrarly large input/output size. Therefore, based on 256 benchmarks we interpolate a formula that returns the performance given an input or a output size. We derive 3 formulas:
- SHA3 hash function targeting 128 bit security. It has variable input size and fixed 256 bit output.
- SHAKE128 XOF with 256 bit input and variable output size.
- SHAKE128 PRG with 128 bit input seed and variable output size.
The three formulas are implemented in the `perf.py` file. They ensure an accurate estimation of the accelerated performance using the corresponding hardware architecture. Unfortunately, these formulas are platform specific, and do cover limited cases (fixed output or input size, fixed target security), yet they are sufficient for estimating signature performance. The point is, if other accelerated benchmarks are to be used with this script, the benchmarking approach is left to implement and deeply depends on the target architecture.


[Saa24] Markku-Juhani O. Saarinen. Accelerating SLH-DSA by two orders of magnitude with a
single hash unit. In Leonid Reyzin and Douglas Stebila, editors, CRYPTO 2024, Part I,
volume 14920 of LNCS, pages 276–304. Springer, Cham, August 2024.