# Rough implementation of Simon's Algorithm using QuTiP

from qutip import *
from itertools import *
from utils import *
import operator
import ualgebra as UA

# TODO Supress warnings about imaginary numbers

# Generates an oracle function for Simon's Algorithm
# In: k, order (?) of the group; D, the hidden subgroup D;
#     c, the common value for f(x) forall x in D
# Out: a function (array-form) that seperates cosets on the hidden subgroup
#     (f(a)=f(b) iff a-b in D)
def gen_oracle(k, D):
    c = 0
    f = [ -1 for _ in range(2**k) ]
    for x in range(2**k):
        for y in range(2**k):
            if x^y in D:
                if x == y and f[x] == -1:
                    f[x] = c
                    c += 1
                else:
                    f[x] = f[y] = f[x] if x < y else f[y]
    return f

def ket_as_list(ket):
    return ket.full().astype(int).flatten().tolist()

def int_to_ket(x, n):
    return tensor([ basis(2, d) for d in int_to_bin(x,n) ])

# Generates an oracle operator for Simon's Algorithm (Note: assumes structure
#     of oracle)
# In: k, the size of the group; f, the oracle function; mult, multiplier for the
#       if applicable
# Out: a (2**n)x(2**n) unitary operator embedding the oracle function
def gen_oracle_op(n, f, arity=2):
    ret = [ [] for _ in range(2**(2*n)) ]
    for x in range(len(f)):
        for y in range(len(f)):
            fx = f[x]
            x_ket = int_to_ket(x,n)
            y_ket = int_to_ket((fx+y)%2**n, n)
            ket = tensor( x_ket, y_ket )
            for i,entry in enumerate(ket_as_list(ket)):
                ret[i].append(entry)
    return ret

# List of registers to preserve (0...n-1)
def ptrace_wrt_regs(obj, ris, n):
    qubits = []
    for i in ris:
        qubits.extend( [i * n + j for j in range(n)] )
    return obj.ptrace(qubits)

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
print(dm)
#dm_to_hist(dm)
