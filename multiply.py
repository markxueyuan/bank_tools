from sympy.solvers.solveset import nonlinsolve
from sympy import Symbol
from sympy import simplify
from sympy.solvers import solve


################## Variables ########################

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
r = Symbol('r')

# Firm

eqF1 = P * r -R_L


# HH, FOC

eq1 = R_D1 * θ * P_B - (R_F0 * θ)
eq2 = R_E1 - R_D1 * θ
eq3 = -2*R_F0*θ + 2*R_E - λ*R_D*θ + R_F0 * θ * λ * R_W / P_B
eq4 = -2*R_F0*θ + R_D*θ + (1-λ)*R_D*θ + R_F0 * θ * λ * R_W / P_B


"""
s1:

[{R_E: R_D*θ, P_B: R_F0*R_W*λ/(R_D*λ - 2*R_D + 2*R_F0)}]

s_1 = solve((eq3, eq4),R_E,P_B, dict=True)

[{R_E: R_D*θ, P_B: R_F0*R_W*λ/(R_D*λ - 2*R_D + 2*R_F0)}]

s_1 = solve((eq3, eq4),R_E, R_F0)

 {R_E: R_D*θ, R_F0: P_B*R_D*(-λ + 2)/(2*P_B - R_W*λ)}

s_1 = solve((eq3, eq4),R_E, R_F0, P_B, dict=True)

[{R_E: R_D*θ, P_B: R_F0*R_W*λ/(R_D*λ - 2*R_D + 2*R_F0)}]

s_1 = solve((eq3, eq4),R_E, R_F0, P_B, R_D, dict=True)
[{R_D: R_E/θ, P_B: R_F0*R_W*θ*λ/(R_E*λ - 2*R_E + 2*R_F0*θ)}]

s_1 = solve((eq3, eq4), R_D, R_E, R_F0, P_B,  dict=True)
[{R_D: R_E/θ, P_B: R_F0*R_W*θ*λ/(R_E*λ - 2*R_E + 2*R_F0*θ)}]

s_1 = solve((eq3, eq4))
[{R_D: R_E/θ, P_B: R_F0*R_W*θ*λ/(R_E*λ - 2*R_E + 2*R_F0*θ)}, {R_E: 0, θ: 0}]

s2:

[{R_E: θ*(R_D*λ - R_D1*R_W*λ + 2*R_F0)/2, R_E1: R_D1*θ, P_B: R_F0/R_D1}]

"""

s_1 = solve((eq3, eq4),R_E, R_F0, λ, R_W, P_B)

s_2 = solve((eq1, eq2, eq3),R_E, R_F0, R_E1, P_B)

s_1 = solve((eq3, eq4))

