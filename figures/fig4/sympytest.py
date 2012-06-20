from sympy import S,solve

s0 = S('s0')
s1 = S('s1')
s2 = S('s2')
kf = S('kf')
kr = S('kr')
L_0 = S('L_0')
R_0 = S('R_0')

eqns = [-s0*s1 + s2, -s0*s1 + s2, s0*s1 - s2, s0 + s2 - 1, s2 + s1 - 1]

eqns2 = [-kf*s0*s1 + kr*s2,
 -kf*s0*s1 + kr*s2,
 kf*s0*s1 - kr*s2,
 -L_0 + s0 + s2,
 -R_0 + s1 + s2]
