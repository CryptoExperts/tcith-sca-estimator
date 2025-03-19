from abc import ABC, abstractmethod
from math import ceil, floor, log2, sqrt, log, pow, comb
from perf import *

def clog2(x) :
    return ceil(log2(x))

# returns the number of keccak permutations given the input/output size and rate (a.k.a. block size)
def keccak_count(input_size, output_size, keccak_rate):
    return ceil((input_size+8)/keccak_rate) + ceil(output_size/keccak_rate) - 1

# return the number of AES calls required to use AES as a PRG in CTR mode given an input/output size and a security parameter
def aes_count(input_size, output_size, security) :
    return ceil(max([input_size, output_size]) / security)

class BinaryTree:
    @staticmethod
    def get_upper_bound_num1(nb_revealed, nb_committed):
        """ Get Upper Bound Num 1:
            Good bound for 'nb_revealed' small """
        return nb_revealed

    @staticmethod
    def get_upper_bound_num2(nb_revealed, nb_committed):
        """ Get Upper Bound Num 2:
            Good bound for 'nb_revealed' big """
        N = nb_committed
        x = nb_revealed
        import math
        return floor((N-x)*math.log2(N/(N-x))) if x < N else 1
        
    @staticmethod
    def get_nb_leaves(nb_revealed, nb_committed, nb_experiments=50):
        """ Get average cost to reveal the seeds """
        import random
        def run_experiment(k, N):
            logN = ceil(log2(N))
            N_ = 2**logN

            arr = [False]*N + [None]*(N_-N)
            for idx in random.sample(list(range(N)), k):
                arr[idx] = True

            count = 0
            for n in [2**i for i in range(logN, 0, -1)]:
                for i in range(0, n, 2):
                    if arr[i+1] is None:
                        arr[i//2] = arr[i]
                        continue
                    if arr[i] != arr[i+1]:
                        count += 1
                    arr[i//2] = arr[i] and arr[i+1]
                for i in range(n//2, n):
                    arr[i] = None
            count += (1 if arr[0] else 0)
            return count

        return sum([
            run_experiment(nb_revealed,nb_committed)
            for _ in range(nb_experiments)
        ]) / nb_experiments

class MPC_protocol :
    def __init__(self, field_size, dim_w, degree=1, kappa=128) :
        self.field_size = field_size
        self.dim_w = dim_w
        self.degree = degree
        self.kappa = kappa
    
    @abstractmethod
    def get_input_size(self):
        pass

    @abstractmethod
    def get_unif_size(self):
        pass

    @abstractmethod
    def get_comm_size(self):
        pass

    @abstractmethod
    def get_randomness_size(self):
        pass
    
    @abstractmethod
    def get_false_positive_rate(self):
        pass

class SDitH_v1_prot(MPC_protocol) :
    def __init__(self, field_size, eta, m, k, w, d, t) :
        super().__init__(field_size, m, 1)
        self.eta = eta  #degree of the field extension F_points
        self.m = m      #code length
        self.k = k      #vector dimension
        self.w = w      #hamming weight bound
        self.d = d      #parameter of the d-splitting variant
        self.t = t      #number of random evaluation points
        self.mpc_rounds = 1

    def get_input_size(self) :
        return (self.k + 2*self.w)*clog2(self.field_size)+self.t*clog2(pow(self.field_size, self.eta))
        
    def get_unif_size(self) :
        return 2*self.d*self.t*clog2(pow(self.field_size, self.eta))

    def get_comm_size(self) :
        return 2*self.d*self.t*clog2(pow(self.field_size, self.eta))

    def get_randomness_size(self) :
        return (1 + self.d)*self.t*self.eta

    def get_false_positive_rate(self):
        # q is the field size
        # m is the code length
        # k is the code dimension
        # w is the weight contraint
        # d is the split factor (d=1 for standard SD instance)
        delta = pow(self.field_size, self.eta)
        split_m = self.m // self.d
        split_w = self.w // self.d

        pr = (split_m+split_w-1)/delta

        # Formula
        p = sum([
            comb(self.t, i) * (pr)**i * (1-pr)**(self.t-i)
            / (delta**(self.t-i))
            for i in range(self.t+1)
        ])
        return p, [p]

class MQOM_v1_prot(MPC_protocol) :
    def __init__(self, field_size, n, n1, n2, eta) :
        super().__init__(field_size, n, 1)
        self.n = n      #number of unknowns
        self.m = n      #number of equations
        self.n1 = n1    #number of coordinates per chunks of x and w.
        self.n2 = n2    #number of chunks in x and w
        self.eta = eta  #extension degree for the field Fqη used in the MPC protocol.
        self.mpc_rounds = 2

    def get_input_size(self) :
        return (self.n + self.eta * (2*self.n1 - 1)) * clog2(self.field_size)
    
    def get_unif_size(self) : 
        return self.n2 * self.eta * clog2(self.field_size)
    
    def get_comm_size(self) :
        return self.n2 * self.eta * clog2(self.field_size)

    def get_randomness_size(self) :
        chal_1 = self.m * self.field_size * self.eta
        chal_2 = self.field_size * self.eta
        return chal_1 + chal_2

    def get_false_positive_rate(self) :
        p1 = 1 / pow(self.field_size, self.eta)
        p2 = (2*self.n1 - 1) / (pow(self.field_size, self.eta) - self.n1)
        p = p1 + (1-p1)*p2
        return p, [p1, p2]

class TCitH_PC_prot(MPC_protocol) :
    def __init__(self, field_size, n, ell=1, rho=None):
        super().__init__(field_size, n, 2)
        self.n = n      #number of unknowns
        self.m = n      #number of equations
        self.ell = ell  #privacy parameter of the TCitH framework
        if rho == None: 
            self.rho = ceil(self.kappa / log2(self.field_size))
        else:
            self.rho = rho
        self.mpc_rounds = 1
        
    def get_input_size(self) :
        return self.n * clog2(self.field_size)

    def get_unif_size(self) :
        d_beta = self.ell
        n_hints = ceil((self.degree-1)*self.ell/(d_beta-self.ell+1))
        return n_hints*self.rho * clog2(self.field_size)
    
    def get_comm_size(self) :
        return self.rho * clog2(self.field_size)

    def get_randomness_size(self) :
        return self.m * clog2(self.field_size)

    def get_false_positive_rate(self) :
        p = 1/pow(self.field_size, self.rho)
        return p, [p]

class MPCitH :
    def __init__(self, N, ell, slack, masking_order, tau=None, mpc_protocol=None, kappa=128, has_aes_prg=False, accel=False, salt_length=None):
        assert kappa in [128, 192, 256]
        self.kappa = kappa
        self.N = N
        self.ell = ell
        self.tau = tau
        self.slack = slack 
        self.nshares = masking_order+1
        self.mpc_protocol = mpc_protocol
        if self.kappa == 128 :
            self.keccak_rate = 1344 # SHAKE128
        elif kappa == 192 or kappa == 256 :
            self.keccak_rate = 1088 # SHAKE256
        self.salt_length = salt_length or self.kappa 
        self.message_length = 0
        self.has_aes_prg = has_aes_prg
        assert self.kappa == 128 or accel == False, 'Hardware acceleration is only supported when kappa=128'
        self.accel=accel

    def get_forgery_cost(self):
        """ Return the forgery cost (in bits) of the KZ20 attack against the scheme """
        N = self.N
        tau = self.tau
        _, p = self.mpc_protocol.get_false_positive_rate()
        t = self.mpc_protocol.mpc_rounds
        if type(p) is not list:
            p = [p]
        assert len(p) == t, (len(p), t)
        ell = self.ell
        if self.with_shamir:
            ell = self.ell
        else:
            ell = 1
        def sum_pmf(tau1, tau, p):
            return sum(
                comb(tau, k)*(pow(p,k))*(pow((1-p),(tau-k)))
                for k in range(tau1, tau+1)
            )
        def rec_cost(tau, p_list, overflow=pow(2,512)):
            if len(p_list) == 1:
                return pow((1/p_list[0]), tau)
            best_cost = None
            for tau1 in range(tau, -1, -1):
                try:
                    c1 = 1 / sum_pmf(tau1, tau, p_list[0])
                except ZeroDivisionError: # Too small value
                    c1 = overflow if (overflow is not None) else pow(2, 512)
                if (overflow is not None) and (c1 >= overflow):
                    # Range 1/3
                    best_cost = overflow
                else:
                    c2 = rec_cost(tau - tau1, p_list[1:])
                    if (overflow is not None) and (c2 >= overflow):
                        # Range 3/3
                        c = c1 + c2
                        if (best_cost is None) or (c < best_cost):
                            best_cost = c
                        break
                    else:
                        # Range 2/3
                        c = c1 + c2
                        if (best_cost is None) or (c < best_cost):
                            best_cost = c
            return best_cost
        # KZ20 Formula
        forgery_cost = rec_cost(tau, p+[1/comb(N,ell)])
        return log(forgery_cost, 2)

    @abstractmethod
    def soundness_error(self):
        pass

    def parallel_repetitions(self) :
        self.tau = 1
        soundness_error = self.soundness_error()
        cost_forge = self.get_forgery_cost()

        acc = soundness_error
        if soundness_error >= 1:
            return -1
            
        goal = pow(2, -self.kappa)
        while(acc > goal or cost_forge < self.kappa) :
            self.tau += 1
            if cost_forge < self.kappa:
                cost_forge = self.get_forgery_cost()
            acc *= soundness_error
        return self.tau

    @abstractmethod
    def get_size(self):
        pass

    @abstractmethod
    def get_keccak_calls(self) :
        pass

    @abstractmethod
    def get_mult_calls_sharing_eval(self) :
        pass

    # return the number of unmasked multiplications used by the computation of the MPC protocol
    def get_unmasked_multiplications(self) :
        total_calls = 0
        if isinstance(self.mpc_protocol, TCitH_PC_prot) :
            total_calls += (self.mpc_protocol.degree*self.ell + 1) * \
            self.mpc_protocol.m * ( (pow(self.mpc_protocol.n, 2)/2) + self.mpc_protocol.n + self.mpc_protocol.rho)
        else: 
            raise NotImplementedError()
        return self.tau * total_calls

    # return the number of masked multiplications used by the computation of the MPC protocol
    def get_masked_multiplications(self) :
        total_calls = 0
        if isinstance(self.mpc_protocol, TCitH_PC_prot) :
            total_calls += (self.mpc_protocol.degree*self.ell + 1) * \
            self.mpc_protocol.m * self.mpc_protocol.n
        else:
            raise NotImplementedError()
        return self.tau * total_calls

    # given an instance, returns the unmasked keccak performance, 
    # depending on the number of calls to keccak for the unaccelerated keccak, 
    # or depending on the input or output size for accelerated keccak 
    def unmasked_cost(self, n_calls_keccak, size_in_out, keccak_perf, accel_perf):
        # when using hw accel, hash/xof performance does not depend on the number of calls to the 
        # permutation but on the input and/or output size of the hash function / xof 
        if self.accel:
            return accel_perf(ceil(size_in_out/8))
        # when not using accel, hash / xof performance depends on the number of calls to the 
        # keccak permutation
        else:
            return n_calls_keccak * keccak_perf[1]