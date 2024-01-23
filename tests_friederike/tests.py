import numpy as np


b=np.array((10,20))
b = np.append(b, np.zeros(3)) if 'b' in locals() else np.array(0)
print(b)
"""

a = np.zeros((1,7))

A = np.vstack([A, a]) if 'A' in locals() else a  # if A exists, add line a, otherwise A=a
A = np.vstack([A, a]) if 'A' in locals() else a  # if A exists, add line a, otherwise A=a

b = np.ones((2,7))
A = np.concatenate((a,b))
print(A)

b=[1,2]
c=[3,4]
for a in b + c:
    print(a)

Ex_C_col = { }
Ex_C_col["therm"] = 3

print(Ex_C_col)

print(type("therm"))



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