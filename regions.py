from sympy.core.symbol import symbols
from sympy.solvers.solveset import nonlinsolve
from sympy import Symbol, simplify
from sympy.solvers import solve
from sympy.solvers.inequalities import solve_rational_inequalities, \
    reduce_inequalities, solve_univariate_inequality

Ba = Symbol("Ba")
Bb = Symbol("Bb")
Ha = Symbol("Ha")
y = Symbol("y")
A = Symbol("A")
B0 = Symbol("B0")
B = Symbol("B")
M = Symbol("M")
W = Symbol("W")
A = Symbol("A")
H0 = Symbol("H0")
I = Symbol("I")

eq1 = Ba + Bb + Ha - y*A
eq2 = B0 - B + M + W - A
eq3 = H0 - W + A

I = y*A - Bb - M

############## Case 8

s_eq1 = eq1.subs(Ba,B0).subs(Bb,B0).subs(Ha,H0)

s8 = solve((s_eq1, eq2, eq3), A, B0, H0, dict=True)
s_B0 = s8[0][B0]
s_H0 = s8[0][H0]
s_A = s8[0][A]
s_Bb = s8[0][B0]
s_I = simplify(I.subs(A,s_A).subs(Bb,s_Bb))

############### Case 10

s_eq1 = eq1.subs(Ba,B0).subs(Bb,B0).subs(Ha,0)

s10 = solve((s_eq1, eq2, eq3), A, B0, H0, dict=True)
s_B0 = s10[0][B0]
s_H0 = s10[0][H0]################# Case 1

s_eq1 = eq1.subs(Ba,0).subs(Bb,0)
s_eq2 = eq2.subs(B0,0)
s1 = solve((s_eq1, s_eq2, eq3), A, H0, Ha, dict=True)
s_H0 = s1[0][H0]
s_Ha = s1[0][Ha]
s_A = s1[0][A]
s_I = simplify(I.subs(A,s_A).subs(Bb,0))
s_A = s10[0][A]
s_Bb = s10[0][B0]
s_I = simplify(I.subs(A,s_A).subs(Bb,s_Bb))

##################### Case 12

s_eq1 = eq1.subs(Ba,0).subs(Bb,B0).subs(Ha,0)

s12 = solve((s_eq1, eq2, eq3), A, B0, H0, dict=True)
s_B0 = s12[0][B0]
s_H0 = s12[0][H0]
s_A = s12[0][A]
s_Bb = s12[0][B0]
s_I = simplify(I.subs(A,s_A).subs(Bb,s_Bb))

################# Case 1

s_eq1 = eq1.subs(Ba,0).subs(Bb,0)
s_eq2 = eq2.subs(B0,0)
s1 = solve((s_eq1, s_eq2, eq3), A, H0, Ha, dict=True)
s_H0 = s1[0][H0]
s_Ha = s1[0][Ha]
s_A = s1[0][A]
s_I = simplify(I.subs(A,s_A).subs(Bb,0))


################# Case 2

s_eq1 = eq1.subs(Ba,0).subs(Bb,0).subs(Ha,H0)
s_eq2 = eq2.subs(B0,0)
s2 = solve((s_eq1, s_eq2, eq3), A, H0, y, dict=True)
s_H0 = s2[0][H0]
s_y = s2[0][y]
s_A = s2[0][A]
s_I = simplify(I.subs(A,s_A).subs(Bb,0))


################# Case 5

s_eq1 = eq1.subs(Bb,B0).subs(Ha,0)
s_eq3 = eq3.subs(H0,0)
s5 = solve((s_eq1, eq2, s_eq3), A, B0, Ba, dict=True)
s_B0 = s5[0][B0]
s_Ba = s5[0][Ba]
s_A = s5[0][A]
s_I = simplify(I.subs(A,s_A).subs(Bb,s_B0))

################# Case 6

s_eq1 = eq1.subs(Ba,B0).subs(Bb,B0).subs(Ha,0)
s_eq3 = eq3.subs(H0,0)
s6 = solve((s_eq1, eq2, s_eq3), A, B0, y, dict=True)
s_B0 = s6[0][B0]
s_y = s6[0][y]
s_A = s6[0][A]
s_I = simplify(I.subs(A,s_A).subs(Bb,s_B0))

