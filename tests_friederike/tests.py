import numpy as np

a=[1,2,3]
print(a[-1])

"""
# don't delete: still important!
A=np.array([[2,1]])
b=np.array([2])
print(A)
print(b)
print(np.linalg.lstsq(A,b)[0])      # for overdetermined matrix use lstsq instead of solve

types = {"therm": 0, "mech": 1, "chemical": 2}
aux_eqs=[[2, "therm"],[3, "mech"],[1, "chemical"]]
aux_eq = aux_eqs[0]
print(types[aux_eq[1]])

a=[0,0,0,1,0]
print(np.any(a))
"""