from tcith_estimator_common import *

class TCitH_GGM(MPCitH) : 
    def __init__(self, N, ell, slack, masking_order, tau=None, mpc_protocol=None, kappa=128, has_d_seed_trees=False, has_aes_prg=False, accel=False) :
        super().__init__(N, ell, slack, masking_order, tau, mpc_protocol, kappa, has_aes_prg, accel)
        self.with_shamir = True
        self.has_d_seed_trees = has_d_seed_trees

        if tau == None:
            self.tau = self.parallel_repetitions()
        if self.tau == -1:
            raise ValueError('Invalid parameter set: soundness error for tau=1 is greater than 1')
        
    def get_size(self) :
        d_alpha = self.mpc_protocol.degree*self.ell
        dig = (self.mpc_protocol.mpc_rounds+1) * 2*self.kappa
        inputs = self.mpc_protocol.get_input_size()
        comm = self.mpc_protocol.get_comm_size()

        if self.has_d_seed_trees:
            commitments = (self.nshares) * self.kappa * clog2(comb(self.N, self.ell - self.slack)) + 2*self.kappa
        else :
            commitments = self.kappa * clog2(comb(self.N, self.ell - self.slack)) + 2*self.kappa

        if isinstance(self.mpc_protocol, TCitH_PC_prot):
            res = dig + self.tau * (inputs + (d_alpha - self.ell + self.slack)*comm + commitments)
            return ceil(res/8)
        else:
            res = dig + self.tau * (inputs + (d_alpha + 1 - self.ell + self.slack)*comm + commitments)
            return ceil(res/8)

    def soundness_error(self) :
        p, _ = self.mpc_protocol.get_false_positive_rate()
        b0 = comb(self.mpc_protocol.degree*self.ell, self.ell-self.slack)
        b1 = comb(self.N, self.ell-self.slack)
        res = p + (1 - p) * (b0/b1)
        return res

    def get_keccak_cost_hash(self, keccak_perf, accel_perf):
        unmasked = 0
        masked = 0

        d_alpha = self.mpc_protocol.degree*self.ell
        size_out = 2*self.kappa

        if self.has_d_seed_trees:
            # number of calls to commit to the masking shares of the seeds individually
            n_calls_hash = comb(self.N, self.ell - self.slack)
            size_in = self.kappa
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.nshares * self.tau * n_calls_hash * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
            # number of calls to compute the commitments
            size_in = self.salt_length + self.nshares * 2*self.kappa
            n_calls_keccak =  keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.tau * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        else:
            # number of calls to commit to the seeds
            n_calls_hash = comb(self.N, self.ell - self.slack)
            size_in = self.kappa
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            if (self.nshares == 1):
                unmasked += self.tau * n_calls_hash * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
            else:
                masked += self.tau * n_calls_hash * n_calls_keccak * keccak_perf[self.nshares]
        # number of calls for the hash of the commitments
        size_in = 2*self.kappa * comb(self.N, self.ell - self.slack) + self.mpc_protocol.get_input_size()
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.tau * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        # number of calls for the first hash h1
        size_in = self.salt_length + self.tau * 2*self.kappa
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        # number of calls for the extra mpc challenges
        for i in range(1, self.mpc_protocol.mpc_rounds) :
            size_in = self.salt_length + (self.tau+1) * 2*self.kappa
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        # number of calls for the last hash h2
        size_in = self.salt_length + self.message_length + 2*self.kappa + self.tau * (d_alpha+1-self.ell)* self.mpc_protocol.get_comm_size()
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        return unmasked, masked 

    def get_keccak_cost_prg(self, keccak_perf, accel_perf):
        unmasked = 0
        masked = 0

        if self.has_d_seed_trees:
            # number of prg calls to build the tree prg 
            n_calls_prg = comb(self.N, self.ell - self.slack)
            size_in = self.kappa
            size_out = 2*self.kappa
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.tau * self.nshares * n_calls_prg * self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
            # number of calls to derive randomness from
            n_calls_prg = comb(self.N, self.ell - self.slack)
            size_in = self.kappa
            size_out = self.mpc_protocol.get_input_size() + self.mpc_protocol.get_unif_size()
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.tau * self.nshares * n_calls_prg * self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
        else:
            # number of prg calls to build the tree prg 
            n_calls_prg = comb(self.N, self.ell - self.slack)
            size_in = self.kappa
            size_out = 2*self.kappa
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            if (self.nshares == 1):
                unmasked += self.tau * n_calls_prg * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
            else:
                masked += self.tau * n_calls_prg * n_calls_keccak * keccak_perf[self.nshares]
            # number of calls to derive randomness from
            n_calls_prg = comb(self.N, self.ell - self.slack)
            size_in = self.kappa
            size_out = self.mpc_protocol.get_input_size() + self.mpc_protocol.get_unif_size()
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            if(self.nshares == 1):
                unmasked += self.tau * n_calls_prg * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
            else:
                masked += self.tau * n_calls_prg * n_calls_keccak * keccak_perf[self.nshares]
        return unmasked, masked

    def get_keccak_cost_xof(self, keccak_perf, accel_perf):
        unmasked = 0
        masked = 0
        
        # number of calls to get randomness for the MPC protocol
        size_in = 2*self.kappa
        size_out = self.tau * self.mpc_protocol.get_randomness_size()                           # to check 
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
        # number of calls to get the set of opened parties
        size_in = 2*self.kappa
        size_out = self.tau * (self.ell - self.slack) * clog2(self.N)
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
        
        return unmasked, masked 


    # compute the required number of multiplications required to compute the parties' secret shares 
    def get_mult_calls_sharing_eval(self):
        umasked = 0
        masked = 0

        unmasked = self.tau * self.ell * ((self.mpc_protocol.get_input_size() + \
                   self.mpc_protocol.get_unif_size()) / clog2(self.mpc_protocol.field_size)) * \
                   (comb(self.N, self.ell-self.slack) + 1)
        unmasked += self.nshares * self.tau * ((self.mpc_protocol.degree*self.ell)+1) * \
                    ((self.mpc_protocol.get_input_size() + \
                    self.mpc_protocol.get_unif_size()) / clog2(self.mpc_protocol.field_size))


        return unmasked, masked

    def get_performance_index(self, keccak_perf, aes_perf, mult_perf, show=False, accel=False) :
        p0_u = 0
        p0_m = 0
        p1 = 0
        p2_u = 0
        p2_m = 0

        unmasked_hash, masked_hash = self.get_keccak_cost_hash(keccak_perf, hash_accel_128)
        unmasked_xof, masked_xof = self.get_keccak_cost_xof(keccak_perf, xof_accel_128)
        unmasked_prg, masked_prg = self.get_keccak_cost_prg(keccak_perf, prg_accel_128)        
        
        # unmasked keccak perf
        p0_u = unmasked_hash + unmasked_prg + unmasked_xof
        # masked keccak perf
        p0_m = masked_hash + masked_prg + masked_xof

        unmasked_mult_pol_eval, _ = self.get_mult_calls_sharing_eval()
        # (unmasked) shares evaluation perf
        p1 += mult_perf[1] * unmasked_mult_pol_eval

        masked = 0
        unmasked = self.get_unmasked_multiplications()
        if ceil((self.nshares) / (self.slack+1)) == 1 : 
            unmasked += self.get_masked_multiplications()
        else:
            masked += self.get_masked_multiplications()

        # MPC unmasked multiplication perf
        p2_u += mult_perf[1] * (ceil((self.nshares) / (self.slack+1))) * unmasked
        # MPC masked multiplication perf
        p2_m += mult_perf[ceil((self.nshares) / (self.slack+1))] * masked

        perf = p0_u + p0_m + p1 + p2_u + p2_m

        detail = tuple([p0_u, p0_m, p1, p2_u, p2_m])
        keccak_detail = tuple([unmasked_hash, masked_hash, unmasked_prg, masked_prg, unmasked_xof])

        return perf, detail, keccak_detail