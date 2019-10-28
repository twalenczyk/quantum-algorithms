# experiment to determine whether the zero congruence class forms the underlying
# set of an abelian group

# imports {{{1
from itertools import *
from sys import *
import ualgebra as UA
import post_ops as PO
import random
#----------------------------------------------------------------------------}}}1

# experiment {{{1
def plus(x,y):
    return [ x[i]*y[i] for i in range(len(x)) ]

Names = [
"T"    , "P0"    , "P1"    , "P"    , "M"   , "MP0"   , "MP1"   , "MP"   ,
"MEET" , "MEETP0", "MEETP1", "MEETP", "JOIN", "JOINP0", "JOINP1", "JOINP",
"D"    , "DP"    , "DM"   , "A"   , "AD"    , "AP0"   , "AP1"   , "AP"   ,
"U"    , "UD"    , "UM"   , "UP0" , "UP1"   , "F"
] + [ "T0k for 2 <= k < inf" ] + [ "T0inf" ]   + \
    [ "PT0k for 2 <= k < inf" ] + [ "PT0inf" ]  + \
    [ "T1k for 2 <= k < inf" ] + [ "T1inf" ]   + \
    [ "PT1k for 2 <= k < inf" ] + [ "PT1inf" ]  + \
    [ "MT0k for 2 <= k < inf" ] + [ "MT0inf" ]  + \
    [ "MPT0k for 2 <= k < inf" ] + [ "MPT0inf" ] + \
    [ "MT1k for 2 <= k < inf" ] + [ "MT1inf" ]  + \
    [ "MPT1k for 2 <= k < inf" ] + [ "MPT1inf" ]


print("Possible clones:")
for name in Names:
    print("  ", name)
inp = input("Input a clone name: ")
print("Clone name", inp, "provided")
Ops = PO.named_clone(inp)
print("Corresponding operations:")
for op in Ops:
    print("  ", op.name) 

n = 3
A = UA.FancySet( initial=[list(a) for a in list(product([0,1],repeat=n))] ) # {0,1}^n

stats = 0
total=10**3
for count in range(1,total+1):
  ng = random.randrange(1, 2**n+1)
  Theta, Gens = UA.rand_cong(A, Ops, num_gen=ng, Progress=False)
  zcong = UA.cong_classes(Theta, A)[0] # Hueristic to improve performance
  is_closed = True
  for (x,y,z) in product(zcong, repeat=3):
    if plus(plus(x,y),z) != plus(x,plus(y,z)):
      is_closed = False
      break
  if is_closed:
    stats += 1
  prog = str(round(stats/count*100, 4))
  stdout.write("\r current progress: " + str(round(count/total*100, 4)) + "%    current success rate: " + prog + "%    ")
  if count % 10 == 0: stdout.flush()
stdout.write("\n success = " + str(round(stats/total*100, 4)) + "%\n")
#----------------------------------------------------------------------------}}}1
