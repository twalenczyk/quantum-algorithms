# Rough implementation of Simon's Algorithm using QuTiP

from qutip import *
from itertools import *
from utils import *
import operator
import ualgebra as UA

# TODO Supress warnings about imaginary numbers




# TODO Implement measurement phase
# Simulates Simon's Algorithm
# In: n, the digits for the binary representation of group elements;
#       U, an oracle operator
# Out: TBD
def SimonsAlg(n,U):
    # Useful structures
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)

    # Prepare the state psi
    psi = ht * zn

    # Apply the oracle
    full_reg = tensor(psi, zn)
    post_oracle = U * full_reg

    # Return to the starting space
    targ = ptrace_wrt_regs(post_oracle, [0], n)
    targ = ht*targ*ht

    return targ

# ~~~ Testing ~~~

def verify_oracle(f, U, n):
    cmp_kets = lambda x,y: ket_as_list(x) == ket_as_list(y)
    for x in range(len(f)):
        for y in range(len(f)):
            # We expect U|xy>=|x,y oplus f(x)>
            fx = f[x]
            yop = (y + fx) % 2**n
            expectation = tensor( int_to_ket(x,n), int_to_ket(yop,n) )

            # Find the reality
            inp = tensor( int_to_ket(x,n), int_to_ket(y,n) )
            reality = U * inp

            if not cmp_kets(expectation, reality):
                print('x y', x,y)
                print('expecation: val ket', yop, expectation)
                print('reality: val ket', '-1', reality)
                print('inp', inp)
                return
    print('Oracle verified')


# Using the partial implementation of Simon's algorithm

n = 3
f = gen_oracle(n, set([0,2]) )

# Programmatically generated operator
op = gen_oracle_op(n,f)
U = Qobj( inpt=op, dims=[[2]*2*n, [2]*2*n])
verify_oracle(f, U, n)

dm = SimonsAlg(n,U)
