# Group implementation using algebra library
# Author: Taylor Walenczyk
# Last updated: 08.30.2019

from utils.summations import dot_prod
from utils.number_representation import *
from utils.printing import *
UA = __import__('ualgebra')
import itertools
import math
import random
import sys
import time

# TODO Refactor to use FancySet class
# TODO Refactor use of elements to accomodate set
class Group(): #{{{
    # Represents groups in ZZ_2**k where k is variable. The group operation at the
    # moment is component wise (i.e. bitwise) multiplication. Future iterations will
    # let this be definable.

    # Functions:
    # __init__(...): Provide the order of the group and the group operation
    #       along with your choice of elements
    # populate_M(self, elements): provide a list of elements with which to populate
    #       the group operation look up table (i.e. M)
    # m(self,x,y): returns the result of the group operation (i.e. M[x][y])

    # TODO Refactor to use FancySet
    def __init__(self, k, gop, e, elements=[]): #{{{
        self.elements = elements
        self.K = k
        self.E = e
        self.gop = gop
        self.Inv = dict()
        self.Op = dict()
        self.populate_Op(self.elements)
    # ========================================================================}}}
    def __len__(self): #{{{
        return len(self.elements)
    # ========================================================================}}}
    def __str__(self): #{{{
        present =   '|G| = {0}\n'.format(len(self))
        present +=  'G = {0}\n'.format(self.elements)
        present +=  'G.E = '+self.e()+'\n'
        present +=  'G.Op = '+pprint_dict(self.Op)
        return present
    # ========================================================================}}}
    # TODO Determine whether or not I need this function
    def populate_Op(self, es): #{{{
        for a in es:
            if a not in self.Op: # Initialize if necessary
                self.Op[a] = dict()
            for b in self.elements:
                if b not in self.Op:
                    self.Op[b] = dict()
                if self.gop == None: # Who knows what to do then?
                    continue
                res = self.gop(a,b)
                add = pad_bin(int_to_bin(res), self.K)
                self.Op[a][b] = add
                self.Op[b][a] = add
                # looking forward to prevent issues
                if add not in self.Op:
                    self.Op[add] = dict()
    # ========================================================================}}}
    def op(self,x,y): #{{{
        return self.Op[x][y]
    # ========================================================================}}}
    def inv(self,x): #{{{
        if x in self.Inv:
            return self.Inv[x]
        for y in self.elements:
            if self.op(x,y) == self.e():
                self.Inv[x] = y
                return y
        return None # if x does not have an inverse
    # ========================================================================}}}
    def e(self): #{{{
        return self.E
    # ========================================================================}}}
    def is_group(self): #{{{
        # idempotence
        for x in self.elements:
            if self.op(x,self.e()) != x:
                return False

        # associativity
        for tup in itertools.product(self.elements,repeat=3):
            x = tup[0]
            y = tup[1]
            z = tup[2]

            try:
                if self.op(self.op(x,y),z) != self.op(x,self.op(y,z)):
                    print('NON-ASSOCIATIVE: x={0}, y={1}, z={2}'.format(x,y,z))
                    return False
            except:
                print('Invalid indexing: x={0}, y={1}, z={2}'.format(x,y,z))
                sys.exit()
        return True
    # ========================================================================}}}
    def is_abelian(self): #{{{
        # commutativity
        for tup in itertools.product(self.elements, repeat=2):
            x = tup[0]
            y = tup[1]

            try:
                if self.Op(x,y) == self.Op(y,x):
                    return False
            except:
                print('Invalid indexing: x={0}, y={1}'.format(x,y))
                sys.exit()
        return False
    # ========================================================================}}}
# ----------------------------------------------------------------------------}}}
# Creates a passable group function
# In:   G, a group
# Out:  a callable function that performs an embedded group's operation
def grp_op(G): #{{{
    def op(x,y):
        return G.op(x,y)
    return op
# ----------------------------------------------------------------------------}}}
# Generate the closure of an operator group
# In:   G, a group; U1 and U2, sets of elements in the operator; z, the added
#       element
# Out:  Nothing, but a correct operator table G.Op
def close_operator(G, U1, U2, ad):
    options = list(U1)
    options.extend(U2)
    s1 = set(U1)

    while len(s1) != 0:
        x = s1.pop()
        a = random.choice(options)
        if ad in G.Op[x]:
            a = G.op(x,ad)
        G.Op[x][ad] = a

        # Form closure of op(op(x,ad),z) = op(a,z) = op(x,b) = op(x,op(ad,z))
        for z in U2:
            b = random.choice(options)
            if z in G.Op[ad]:
                b = G.op(ad,z)
            G.Op[ad][z] = b

            c = random.choice(options)
            if z in G.Op[a] and b in G.Op[x]:
                if G.op(a,z) != G.op(x,b):
                    print('ERROR: Non-associativity detected when closing operator')
                    print('\t m.op({0},{1})={2} while m.op({3},{4})={5}'.format(
                        a,z,G.op(a,z),x,b,G.op(x,b)))
                    sys.exit()
                c = G.op(a,z)
            elif z in G.Op[a]:
                c = G.op(a,z)
            elif b in G.Op[x]:
                c = G.op(x,b)

            G.Op[a][z] = c
            G.Op[x][b] = c
# ----------------------------------------------------------------------------}}}
# Generates a random group
# In:   Some stuff
# Out:  A random group
def rand_group(): #{{{
    # Pick a random number of digits to use
    digs = random.randrange(1,3+1)
    elements = [ ''.join(tup) for tup in itertools.product(['0','1'], repeat=digs) ]
    g = Group(
                k=digs,
                gop=None,
                e=elements[ random.randrange(len(elements)) ],
                elements=elements
            )
    # Fill out identity column
    e = g.e()
    for x in g.elements:
        g.Op[e][x] = x
        g.Op[x][e] = x
    # Randomly fill in multiplication table
    # Pick two elements, randomly assign their operation
    s1 = set(g.elements)
    u1 = []
    u2 = []
    while len(s1) != 0:
        x = s1.pop()
        u1.append(x)
        s2 = set(g.elements)
        while len(s2) != 0:
            y = s2.pop()
            u2.append(y)
            a = random.choice(g.elements)
            if y in g.Op[x]:
                a = g.op(x,y)
            g.Op[x][y] = a

            # Form closure of op(op(x,y),z) = op(a,z) = op(x,b) = op(x,op(y,z))
            for z in g.elements:
                b = random.choice(g.elements)
                if z in g.Op[y]:
                    b = g.op(y,z)
                g.Op[y][z] = b

                c = random.choice(g.elements)
                if z in g.Op[a] and b in g.Op[x]:
                    if g.op(a,z) != g.op(x,b):
                        print('ERROR: Non-associativity detected when generating random group.')
                        print('\t m.op({0},{1})={2} while m.op({3},{4})={5}'.format(
                            a,z,g.op(a,z),x,b,g.op(x,b)))
                        sys.exit()
                    c = g.op(a,z)
                elif z in g.Op[a]:
                    c = g.op(a,z)
                elif b in g.Op[x]:
                    c = g.op(x,b)

                g.Op[a][z] = c
                g.Op[x][b] = c
                close_operator(g, u1, u2, c)

    return g
# ----------------------------------------------------------------------------}}}
# Generates a proper group from generating elements
# In:   G, a group; gens, the generating elements
# Out:  a FancySet of elem ents in the subgroup
def gen_subgroup(G, gens): #{{{
    cur_set = comp_lib.FancySet(initial=gens, addl=[ 'generator' for _ in gens ])
    new_set = comp_lib.FancySet(initial=[G.e()], addl='identity element') # ensures e exists
    # Add inverse elements
    for x in cur_set:
        xi = G.inv(x)
        if xi and xi not in cur_set: # prevents overwriting addl
            new_set.add(xi, addl='Inverse of {0}'.format(x))

    while len(new_set) != 0: # while the closure grows
        cur_set.update(new_set)
        new_set = comp_lib.FancySet()
        for ps in itertools.product(cur_set, repeat=2):
            res = G.op(ps[0],ps[1])
            if res not in cur_set:
                note = 'Generated by G.op({0},{1})'.format(ps[0],ps[1])
                new_set.add(res, addl=note)
        # MAY NEED TO WRAP THIS IN ITS OWN CLOSURE
        # Add inverse elements
        for x in new_set: # Presumes inverse are present already for old set
            xi = G.inv(x)
            if xi and xi not in new_set: # prevents overwriting addl
                new_set.add(xi, addl='Inverse of {0}'.format(x))
        # Add elements to ensure assoc
        new_union = cur_set.union(new_set)
        for tup in itertools.product(new_union, repeat=3):
            x = tup[0]
            y = tup[1]
            z = tup[2]

            a = G.op(x,y)
            b = G.op(y,z)
            c = G.op(a,z)

            if a not in cur_set:
                new_set.add(a)
            if b not in cur_set:
                new_set.add(b)
            if c not in cur_set:
                new_set.add(c)
        # Add inverse elements
        for x in new_set:
            xi = G.inv(x)
            if xi and xi not in new_set: # prevents overwriting addl
                new_set.add(xi, addl='Inverse of {0}'.format(x))
    return cur_set
# ----------------------------------------------------------------------------}}}
# Creates a random subgroup
# In:   num_gens, the number of generators for the subgroup (-1 if a random amount
#       is desired); G, a group object; Ops, list of relevant operators
# Out:  a subgroup of G
def rand_subgroup(G, num_gens=-1): #{{{
    if num_gens == -1:
        num_gens = random.randrange(len(G))
    gens = comp_lib.FancySet()
    for _ in range(num_gens):
        gens.add(random.choice(G.elements), addl='generator')
    S = gen_subgroup(G, gens)
    return S
# ----------------------------------------------------------------------------}}}
# =============================== Tests =========================================
def test_rand_subgroups(iterations=1000, digits=5): #{{{
    iters = iterations
    digs = digits
    G = Group(
                k=digs,
                gop=lambda x,y: int(x,2) & int(y,2),
                e='1'*digs,
                elements=[ ''.join(tup) for tup in itertools.product( ['0', '1'], repeat=digs ) ]
            )
    for _ in range(iters):
        S = rand_subgroup(G)
        gs = Group(
                    k=G.K,
                    gop=None,
                    e=G.e(),
                    elements= [ i for i in S ] # Snags the elements
                )
        gs.Op = dict()
        # Restrict the operator table to ensure associativity
        for x in S:
            gs.Op[x] = G.Op[x]
            for y in G.Op[x]:
                if y not in G.Op[x]:
                    del gs.Op[x][y]
        if not gs.is_group():
            print('================ FAILURE ===============')
            print(S)
            print(gs)
            print('================ FAILURE ===============')
            return
# ----------------------------------------------------------------------------}}}
if __name__ == '__main__': #{{{
    g = Group(
            k=5,
            gop=lambda x,y: int(x,2) & int(y,2),
            e='1'*5,
            elements=[ ''.join(tup) for tup in itertools.product( ['0','1'], repeat=5 ) ]
        )
    bitwise_mult = grp_op(g)
    M = UA.Operation(bitwise_mult, 2, 'bitwise multiplication')

    start_time = time.time()
    #test_rand_subgroups()
    rsg = UA.rand_subgrp(g, [M])
    end_time = time.time()
    print(str(rsg))
    print('Is this a group? {0}'.format(rsg.is_group()))
    print('Total execution time: {0} seconds'.format(end_time-start_time))
# ----------------------------------------------------------------------------}}}
