from sympy.solvers.solveset import nonlinsolve
from sympy import Symbol, symbols, init_printing,\
    simplify, integrate, latex, solve_poly_system,\
    Function, Derivative, apart, cancel, expand, collect,\
    Integral, pprint, srepr

from sympy.solvers import solve, nonlinsolve

import pdb


# generate prettier output
init_printing()


θ = Symbol('θ')
R_E1 = Symbol('R_E1')
R_D1 = Symbol('R_D1')
μ_D1 = Symbol('μ_D1')
μ_E1 = Symbol('μ_E1')
μ_Bh = Symbol('μ_Bh')
R_B = Symbol('R_B')
P_B = Symbol('P_B')
R_E = Symbol('R_E')
R_D = Symbol('R_D')
R_W = Symbol('R_W')
R_F0 = Symbol('R_F0')
μ_λa_ = Symbol('μ_λa_')
μ_λb_ = Symbol('μ_λb_')
μ_λa = Symbol('μ_λa')
μ_λb = Symbol('μ_λb')
D = Symbol('D')
E = Symbol('E')
A = Symbol('A')
W = Symbol('W')
L = Symbol('L')
M = Symbol('M')
D1 = Symbol('D1')
E1 = Symbol('E1')
A1 = Symbol('A1')
B_ha = Symbol('B_ha')
B_hb = Symbol('B_hb')
Bh = Symbol('Bh')
B_B = Symbol("B_B")
λa = Symbol('λa')
λb = Symbol('λb')
λ = Symbol('λ')
α = Symbol('α')
Ma = Symbol('Ma')
Mb = Symbol('Mb')
R_M = Symbol('R_M')
Ba = Symbol('Ba')
Bb = Symbol('Bb')
I = Symbol('I')
Pf = Symbol('Pf')
R_I = Symbol('R_I')
R_L = Symbol('R_L')
μ_ma = Symbol('μ_ma')
μ_Ba = Symbol('μ_Ba')
μ_Bb = Symbol('μ_Bb')
μ_B = Symbol('μ_B')
P_B1 = Symbol('P_B1')
B = Symbol('B')
P = Symbol('P')
e0 = Symbol('e0')
e1 = Symbol('e1')
F = Symbol('F')
B_H = Symbol('B_H')
B1 = Symbol('B1')


############# Equations

eq1 = R_D1 * θ * P_B - (R_F0 * θ)
eq2 = R_E1 - R_D1 * θ
eq4 = -2*R_F0*θ + R_D*θ + (1-λ)*R_D*θ + R_F0 * θ * λ * R_W / P_B

s_e0 = 2*α*A
s_e1 = (1+2*λ*R_W)*α*A

beq3 = R_M - R_D1*(1-(α*(A+A1)+A1*α)/R_E1)-R_E1*(α*(A+A1)+A1*α)/R_E1
beq4 = 2*R_L - 2*R_E*s_e0/R_E - R_D*(1-s_e0/R_E)- R_D1*(A1*α/R_E1) \
       - R_E1*(A1*α/R_E1)



#######################################################
############ 1. replicate Wenting's result ############
#######################################################

s2_R_E = solve(eq2,R_E1)[0]

beq31 = beq3.subs(R_E1, s2_R_E)

s_R_D1 = solve(beq31,R_D1)[0]

beq41 = beq4.subs(R_E1, s2_R_E)

s_R_D = solve(beq41,R_D)[0]

s2_R_F0 = solve(eq1, R_F0)[0]

eq41 = eq4.subs(R_F0, s2_R_F0)

eq42 = eq41.subs(R_D,s_R_D).subs(R_D1,s_R_D1)

s_P_B = solve(eq42, P_B)[0]

s_R_F0 = s2_R_F0.subs(P_B,s_P_B).subs(R_D1,s_R_D1)

s_R_E1 = s2_R_E.subs(R_D1,s_R_D1)


for i in ["s_R_D","s_R_D1","s_R_E1","s_P_B","s_R_F0"]:
    print(i, ": ", eval(i).free_symbols)



"""
The above procedures solves the variables as expressions
of variables in the brackets.

s_R_D :  {A1, A, θ, R_L, α, R_E}
s_R_D1 :  {A1, A, α, R_M, θ}
s_R_E1 :  {A1, A, α, R_M, θ}
s_P_B :  {A1, A, θ, R_M, R_L, R_W, α, R_E, λ}
s_R_F0 :  {A1, A, θ, R_M, R_L, R_W, α, R_E, λ}

I adopt a method slightly different from Wenting's:

Instead of partially substitute R_F0 in eq4, I subsitute all

R_F0 to P_B*R_D1 and solve P_B before solving R_F0.

"""

##############################################################
############ 2. stepping through solve() function ############
##############################################################



def tt():
    pdb.set_trace()
    solve((eq1, eq2, eq4, beq3, beq4),R_D, R_D1, P_B, R_F0, R_E1)


"""

I use python's pdb debugging tool to trace every step of solve()

calculation.

Finally I find that it uses Groebner basis approach:

    Returns all possible solutions over C[x_1, x_2, ..., x_m] of a
    set F = { f_1, f_2, ..., f_n } of polynomial equations,  using
    Groebner basis approach. For now only zero-dimensional systems
    are supported, which means F can have at most a finite number
    of solutions.

    The algorithm works by the fact that, supposing G is the basis
    of F with respect to an elimination order  (here lexicographic
    order is used), G and F generate the same ideal, they have the
    same set of solutions. By the elimination property,  if G is a
    reduced, zero-dimensional Groebner basis, then there exists an
    univariate polynomial in G (in its last variable). This can be
    solved by computing its roots. Substituting all computed roots
    for the last (eliminated) variable in other elements of G, new
    polynomial system is generated. Applying the above procedure
    recursively, a finite number of solutions can be found.

    The ability of finding all solutions by this procedure depends
    on the root finding algorithms. If no solutions were found, it
    means only that roots() failed, but the system is solvable. To
    overcome this difficulty use numerical algorithms instead.

In our problem, the program gets 7 Groebner basis, but the root
finding algorithm failed. But it doesn't mean
the system is not solvable.

For details of the implementation, see:

http://docs.sympy.org/dev/_modules/sympy/solvers/polysys.html#solve_generic


Take out R_FO and take in R_E into our 5X5 problem:
solve((eq1, eq2, eq4, beq3, beq4),R_D, R_D1, P_B, R_E1, R_E, dict=True)
The problem solves.

Groebner basis approach is mathematically sophisticated, we cannot
fully solve the problem unless we can understand the maths.

"""




##############################################################
############ 3. Number of variables               ############
##############################################################


"""
1. When we pass no variables to solve functions, what it does is

to take the union of free symbols of all equations as variables.

2. Groebner algorithm itself doesn't requires the number of
variables should be equal or less than the number of equations

3. The rule of discarding excess solutions:
If depends on previously solved variables: discard it.

I couldn't follow the details of its implementation,
the implementation is on line 1708 of the below page:
http://docs.sympy.org/dev/_modules/sympy/solvers/solvers.html#_solve_system


"""




