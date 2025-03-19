from tcith_estimator_common import *

class TCitH_MT(MPCitH) :
    def __init__(self, N, ell, slack, masking_order, tau=None, mpc_protocol=None, kappa=128, eta=None, has_pr_shares=False, has_pr_shares_1seed=False, has_aes_prg=False, accel=False) :
        super().__init__(N, ell, slack, masking_order, tau, mpc_protocol, kappa, has_aes_prg, accel)
        assert mpc_protocol != None
        assert N <= mpc_protocol.field_size
        if eta == None:
            self.eta = ceil(log(comb(self.N, self.ell+2)*pow(2, self.kappa), self.mpc_protocol.field_size))
        else:
            self.eta = eta
        self.with_shamir = True
        self.has_pr_shares = has_pr_shares
        self.has_pr_shares_1seed = has_pr_shares_1seed
        if tau == None:
            self.tau = self.parallel_repetitions()
        if self.tau == -1:
            raise ValueError('Invalid parameter set: soundness error for tau=1 is greater than 1')
        
    def get_size(self) :
        d_alpha = self.mpc_protocol.degree*self.ell
        d_beta = self.ell

        dig = (self.mpc_protocol.mpc_rounds+2) * 2*self.kappa
        inputs = self.mpc_protocol.get_input_size()
        unif = self.mpc_protocol.get_unif_size()
        comm = self.mpc_protocol.get_comm_size()
        commitments = 2*self.kappa * (self.ell - self.slack) * clog2(self.N / (self.ell - self.slack))
        #commitments = 2*self.kappa * BinaryTree.get_nb_leaves(self.N-(self.ell-self.slack),self.N)
        deg_enforcing = (d_beta+1)*self.eta*clog2(self.mpc_protocol.field_size)
        if self.has_pr_shares:
            pr_seeds = (self.nshares-1)*self.kappa*(self.ell - self.slack)
        elif self.has_pr_shares_1seed :
            pr_seeds = self.kappa*(self.ell - self.slack)
        else :
            pr_seeds = 0

        if isinstance(self.mpc_protocol, TCitH_PC_prot):
            res = dig + self.tau * (
            (self.ell - self.slack) * (inputs + unif) +
            (d_alpha - self.ell + self.slack) * comm +  
            commitments + deg_enforcing + pr_seeds)
            return ceil(res/8)
        else:
            res = dig + self.tau * (
            (self.ell - self.slack) * (inputs + unif) +
            (d_alpha + 1 - self.ell + self.slack) * comm +  
            commitments + deg_enforcing + pr_seeds)
            return ceil(res/8)

    def soundness_error(self) :
        p, _ = self.mpc_protocol.get_false_positive_rate()
        b0 = comb(self.mpc_protocol.degree*self.ell, self.ell-self.slack)
        b1 = comb(self.N, self.ell-self.slack)
        b2 = comb(self.N, self.ell + 2) / pow(self.mpc_protocol.field_size, self.eta)
        res = p + (1 - p) * (b0/b1) + b2
        return res

    def get_keccak_cost_hash(self, keccak_perf, accel_perf): 
        d_alpha = self.mpc_protocol.degree*self.ell
        d_beta = self.ell

        unmasked = 0
        masked = 0
        size_out = 2*self.kappa

        # COMMITMENTS, NO TWEAK => MASKED, EXCEPT WITH APPROPRIATE SLACK VALUE
        if (not self.has_pr_shares) and (not self.has_pr_shares_1seed) :
            # number of calls to commit to each party's share
            size_in = self.mpc_protocol.get_input_size() + self.mpc_protocol.get_unif_size() + clog2(self.mpc_protocol.field_size) * self.eta
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            if (ceil((self.nshares) / (self.slack+1))-1) == 0: 
                unmasked += self.tau * self.N * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
                #print("commitments", self.tau * self.N * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf))
            else:
                masked += self.tau * self.N * n_calls_keccak * keccak_perf[ceil((self.nshares) / (self.slack+1))]
        # COMMITMENTS, TWEAKS => UNMASKED
        else: 
            # number of calls to tompute the intermediate commitments for the first masking share
            size_in = (self.mpc_protocol.get_input_size() + self.mpc_protocol.get_unif_size() + clog2(self.mpc_protocol.field_size) * self.eta)     
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.tau * self.N * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
            # number of calls to tompute the intermediate commitments for the others masking shares
            size_in = self.kappa    
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += (ceil((self.nshares) / (self.slack+1))-1) * self.tau * self.N * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
            # number of calls to compute the output commitments per party 
            size_in = self.salt_length + (ceil((self.nshares) / (self.slack+1))) * 2*self.kappa
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.tau * self.N * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        # number of calls to compute the Merkle Tree
        size_in = 4*self.kappa
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.tau * (self.N-1) * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        #print("merkle tree", self.tau * (self.N-1) * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf))
        # number of calls for the first hash h1
        size_in = self.salt_length + self.tau * 2*self.kappa
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        #print("h1", self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf))
        # number of calls for the degree-enforcing commitment hash 
        size_in = d_beta * self.eta * clog2(self.mpc_protocol.field_size)                           # to check 
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.tau * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        #print("dec hash", self.tau * self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf))
        # number of calls for the second hash h1'
        size_in = self.salt_length + (self.tau+1) * 2*self.kappa
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        #print("h1'",self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf))
        # number of calls for the extra mpc challenges
        for i in range(1, self.mpc_protocol.mpc_rounds) :
            size_in = self.salt_length + (self.tau+1) * 2*self.kappa
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        # number of calls for the last challenge h2
        size_in = self.salt_length + self.message_length + 2*self.kappa + self.tau * (d_alpha+1-self.ell)* self.mpc_protocol.get_comm_size()
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf)
        #print("h2",  self.unmasked_cost(n_calls_keccak, size_in, keccak_perf, accel_perf))
        #print("hash", unmasked, masked)
        return unmasked, masked

    def get_keccak_cost_xof(self, keccak_perf, accel_perf): 
        unmasked = 0
        masked = 0

        size_in = 2*self.kappa
        # number of calls to sample matrices for the degree enforcing commitment scheme
        size_out = self.tau * (self.eta * self.mpc_protocol.dim_w * clog2(self.mpc_protocol.field_size))
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
        #print("matrix", self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf))
        # number of calls to get randomness for the MPC protocol
        size_out = self.tau * self.mpc_protocol.get_randomness_size()                                       # to check 
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
        #print("MPC", self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf))
        # number of calls to expand the opened parties
        size_out = self.tau * (self.ell - self.slack) * clog2(self.N)
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
        #print("parties", self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf))
        return unmasked, masked

    def get_keccak_cost_prg(self, keccak_perf, accel_perf):
        unmasked = 0
        masked = 0

        d_beta = self.ell
        # number of calls to derive randomness from
        size_in = self.salt_length + self.kappa
        size_out = (d_beta+1) * ((self.mpc_protocol.get_input_size() + self.mpc_protocol.get_unif_size()) + self.eta * clog2(self.mpc_protocol.field_size))
        n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
        unmasked += self.tau * self.nshares * self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
        #print("randomness", self.tau * self.nshares * self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf))
        # number of calls for the masks compression tweak
        if self.has_pr_shares or self.has_pr_shares_1seed:
            size_in = self.kappa
            size_out = (self.mpc_protocol.get_input_size() + self.mpc_protocol.get_unif_size()) \
                     + self.eta * clog2(self.mpc_protocol.field_size)
            n_calls_keccak = keccak_count(size_in, size_out, self.keccak_rate)
            unmasked += self.tau * (ceil((self.nshares) / (self.slack+1))-1) * self.N * self.unmasked_cost(n_calls_keccak, size_out, keccak_perf, accel_perf)
        return unmasked, masked

    # compute the required number of multiplications required to compute the parties' secret shares 
    def get_mult_calls_sharing_eval(self):
        unmasked = 0
        masked = 0
        unmasked = self.tau * self.N * self.ell * (((self.mpc_protocol.get_input_size() + \
                   self.mpc_protocol.get_unif_size()) + self.eta * \
                   clog2(self.mpc_protocol.field_size)) / clog2(self.mpc_protocol.field_size))
        return unmasked, masked

    
    # unused currently
    def get_masked_prg_calls_tweak(self):
        nshares = ceil((self.nshares) / (self.slack+1))
        calls = {}
        for i in range(0, 33):
            calls[i] = 0
        def _rec(k) :
            if k == 1 :
                return 
            else :
                c = ceil(k/2)
                f = floor(k/2) 
                calls[k] += self.tau * self.N * keccak_count(self.kappa, 2*self.kappa, self.keccak_rate)
                _rec(c)
                _rec(f)
        if self.has_pr_shares_1seed :
            _rec(nshares)
        return calls
    
    # compute the number of word_len bits refresh calls required by the pseudorandom shares tweak
    def get_refresh_calls(self):
        # |w| + |beta| + |u|
        # word size in bits 
        word_len = 64
        masked = 0
        if self.has_pr_shares or self.has_pr_shares_1seed :
            masked += self.tau * self.N * ceil((self.eta*self.mpc_protocol.field_size + \
                      self.mpc_protocol.get_input_size() + self.mpc_protocol.get_unif_size()) / word_len)
        return masked

    def get_performance_index(self, keccak_perf, aes_perf, mult_perf, refresh_perf, show=False, accel=False) :
        p0_u = 0
        p0_m = 0
        p1 = 0
        p2_u = 0
        p2_m = 0
        p3 = 0
        p4 = 0
        
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

        calls = self.get_masked_prg_calls_tweak()
        
        # not used currently 
        for i in keccak_perf :
            p3 += keccak_perf[i] * calls[i]
            
        masked = self.get_refresh_calls()
        # refresh gadget perf
        p4 += refresh_perf[ceil((self.nshares) / (self.slack+1))] * masked

        perf = p0_u + p0_m + p1 + p2_u + p2_m  + p3 + p4

        keccak_detail = tuple([unmasked_hash, masked_hash, unmasked_prg, unmasked_xof])
        detail = tuple([p0_u, p0_m, p1, p2_u, p2_m, p3, p4])
        return perf, detail, keccak_detail