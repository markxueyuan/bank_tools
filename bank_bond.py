from sympy.core.symbol import symbols
from sympy.solvers.solveset import nonlinsolve
from sympy import Symbol
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
eq1 = R_E1 / θ - R_D1 - μ_D1 + μ_E1
eq2 = R_B - P_B * R_D1 * (1 + μ_D1)
eq3 = R_E / θ - R_D - .5 * (μ_λa_ + μ_λb_)
eq4 = R_D - R_W * (R_D1 + μ_D1) - μ_λa / (D + E) + μ_λa_
eq5 = R_D - R_W * R_D1 *(1 + μ_D1) - μ_λb / (D + E) + μ_λb_
eq6 = (2- λa- λb)-2*R_B + R_W * R_D1* (λa + λb*(1+μ_D1)) \
    + μ_D1*R_W*λa - μ_Bh + μ_λa_*(1-λa) + μ_λb_*(1-λb)

# HH, Slackness condition
eq7 = μ_E1 * E1
eq8 = μ_Bh * (W - (D+E))
eq9 = μ_λa * λa
eq10 = μ_λb * (λb - λ)
eq11 = μ_D1 * (λa * (D + E) * R_W + B_ha*P_B-E1)
eq12 = μ_λa_*(D-λa*(D+E))
eq13 = μ_λb_*(D-λb*(D+E))

# HH, Definitions
eq14 = Bh - W + D + E
eq15 = B_hb * P_B + λb * (D + E)*R_W
eq16 = B_ha * P_B -(D1+E1) + λa*(D+E)*R_W

# Bank, Definitions

eqB1 = D - A + E
eqB2 = D1 - A1 + E1
eqB3 = L - A + M + B_B
eqB4 = E - α * A * A / R_E
eqB5 = E1 - α * A1 * (A + A1) / R_E1
eqB6 = Ma - R_M * M - Ba * P_B + I + λa * A * R_W - A1
eqB7 = Mb - R_M * M - Bb * P_B - I + λa * A * R_W

# Bank, FOC

eqB8 = R_I - R_M - Pf * I - μ_ma
eqB9 = 2 * R_L - R_M*R_M - R_M * R_I - μ_ma * R_M
eqB10 = R_L - R_B - μ_B - μ_Ba - μ_Bb
eqB11 = R_B - P_B * R_M - μ_ma * P_B + μ_Ba
eqB12 = R_B - P_B * R_I + μ_Bb
eqB13 = R_D1 - R_M + (θ-1)/θ * α * (A + 2*A1) - μ_ma
eqB14 = (2-λa-λb)*R_D - 2*R_L + R_W*(λa*R_M + λb*R_I) \
        + (θ-1)/θ * α * (4*A + A1) + μ_ma * λa * R_W


# Bank, slackness conditions

eqB15 = μ_ma * (R_M * M + Ba * P_B1 + A1 - I - λa*A*R_W)
eqB16 = μ_B * B
eqB17 = μ_Ba * (B_B - Ba)
eqB18 = μ_Bb * (B_B - Bb)


# market clearing conditions

eqC1 = B_B + Bh + M - .5* B
eqC2 = Ba + Bb + B_ha + B_hb
eqC3 = Ma + Mb - 2*R_M*M
eqC4 = I - (λb*A*R_W - R_M * M - Bb * P_B)



# redundant variables

eqR1 = P_B1 - P_B









s_HH = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9
            ,eq10, eq11, eq12, eq13)
           ,R_E1, R_D1, μ_D1, μ_E1, μ_Bh, R_B, P_B, R_E,
           R_D, μ_λa_, μ_λb_, μ_λa, μ_λb, B_ha)


s_HH = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9
            ,eq10, eq11, eq12, eq13)
           ,R_E1, R_D1, μ_E1, μ_Bh, R_B, P_B, R_E,
           R_D, μ_λa_, μ_λb_, μ_λa, μ_λb, B_ha,
             dict=True)



s_Bank = solve((eqB1, eqB2, eqB3, eqB4, eqB5, eqB6, eqB7,
            eqB8, eqB9, eqB10, eqB11, eqB12, eqB13,
            eqB14, eqB15, eqB16, eqB17, eqB18),
               R_I, R_L, R_D1, R_D,
               μ_ma, μ_Ba, μ_Bb, μ_B)



s4 = solve((eqF1,
            eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9,
            eq10, eq11, eq12, eq13, eq14, eq15, eq16,
            eqB1, eqB2, eqB3, eqB4, eqB5, eqB6, eqB7,
            eqB8, eqB9, eqB10, eqB11, eqB12, eqB13,
            eqB14, eqB15, eqB16, eqB17, eqB18,
            eqC1, eqC2, eqC3, eqC4,
            eqR1),
           R_E1, R_D1, R_B, P_B, R_E, R_D, R_W, R_I, R_L, P_B1, P,
           D, E, A, L, D1, E1, A1, I,
           B_ha, B_hb, Bh, B_B, Ba, Bb,
           λa, λb, Ma, Mb,
           μ_D1, μ_E1, μ_Bh, μ_λa_, μ_λb_, μ_λa, μ_λb,
           μ_ma, μ_Ba, μ_Bb, μ_B)








































