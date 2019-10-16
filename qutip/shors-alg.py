# Rough implementation of the period finding portion
# of Shor's factoring Algorithm using QuTiP

from qutip import *
from itertools import *
from utils import *
from random import *
import operator
import ualgebra as UA
import math
from cmath import *

# TODO Supress warnings about imaginary numbers

def gen_ua_oracle(n,q,a):
    # Recall that the function acts trivially on vectors >= q
    f = [ x for x in range(2**n) ]
    for x in range(q):
        f[x] = (a*x) % q
    return f

def gen_qft_op(n):
    structure = [ [1] * 2**n ] * 2**n
    exp = lambda a,b : e**(2*pi*a*b*1j/(2**n))
    #complex(0,1) === 0+1j
    for a in range(2**n):
        for b in range(2**n):
            structure[a][b] = exp(a,b)
    return Qobj(
                inpt=structure,
                dims=[ [2]*n, [2]*n ]
            )

def iqft_circuit(n, state):
    # swap qubits
    s = Qobj(
                inpt = state.dims[::-1],
                dims = state.dims
            )
    iqft = gen_qft_op.trans()
    s = iqft*s
    return s

def period_finding(n,q,a):
    f = gen_ua_oracle(n,q,a)
    U = gen_oracle_op(n,f)

    # perform the algorithm!
    zn = tensor([ basis(2, 0) for _ in range(n) ])
    ht = hadamard_transform(n)
    xi = ht * zn
    fxi = U*xi

    # Perform inverse quantum Fourier Transform
    post_qft = iqft_circuit(n, fxi)

    # Discard unnecessary information
    dm = ptrace_wrt_regs(post_qft,[0],n)

    return dm

def simons_alg(y):
    # TODO Change this
    q = randrange(1,y)

    # Check for parity
    if y % 2 == 0:
        return 2

    # Check whether y is the k-th power of an integer for k=[2...lg(y)].
    # If y = m^k, then give the answer m.
    # TODO make better
    for k in range(2, math.ceil(math.log2(y))+1):
        a = math.log(y,k)
        ap = math.floor(a)
        if a-ap == 0: # integer multiple of m^k
            return int(a)

    # Choose an integer a randomly and uniformly between 1 and y-1
    # Compute b=gcd(a,q) (say, by Euclid's algorithm). If b>1, then
    # give answer b.
    a = randrange(1, y)
    b = math.gcd(a,q)
    if b > 1:
        return b

    # Run the quntum algorithm
    r = 1
    if r % 2 == 1:
        return -1

    # Final step
    d = gcd(a**(int(r/2))-1, y)
    if d > 1:
        return d

    return -1


# ~~~ Testing ~~~

# Using the partial implementation of Simon's algorithm

n = 2
Dn = [list(a) for a in list(product([0,1],repeat=n))]

y = randrange(1, 2**n)

op = gen_qft_op(n)
for row in op:
    print(row)



#p1 = simons_alg(y)

#print("(y,p1)=>",(y,p1))

#f = gen_ua_oracle(n,)
#print(f)

#structure = gen_oracle_op(n, f)
#Phi = Qobj( inpt=structure, dims=[[2]*2*n, [2]*2*n] )
#print(Phi)

#dm = clone_alg(n,Phi)
#print("DM trace", dm.tr())
#dm_to_hist(dm)
