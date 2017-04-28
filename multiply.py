from sympy.solvers.solveset import nonlinsolve
from sympy import Symbol, symbols, init_printing,\
    simplify, integrate, latex, solve_poly_system,\
    Function, Derivative, apart, cancel, expand, collect,\
    Integral, pprint

from sympy.solvers import solve, nonlinsolve


# generate prettier output
init_printing()


################## Variables ########################

# symbols defined with assumptions

θ = Symbol('θ', real=True, positive=True)
R_E1 = Symbol('R_E1', real=True)
R_D1 = Symbol('R_D1', real=True)
μ_D1 = Symbol('μ_D1', real=True)
μ_E1 = Symbol('μ_E1', real=True)
μ_Bh = Symbol('μ_Bh', real=True)
R_B = Symbol('R_B', real=True)
P_B = Symbol('P_B', real=True)
R_E = Symbol('R_E', real=True)
R_D = Symbol('R_D', real=True)
R_W = Symbol('R_W', real=True)
R_F0 = Symbol('R_F0', real=True)
μ_λa_ = Symbol('μ_λa_', real=True)
μ_λb_ = Symbol('μ_λb_', real=True)
μ_λa = Symbol('μ_λa', real=True)
μ_λb = Symbol('μ_λb', real=True)
D = Symbol('D', real=True)
E = Symbol('E', real=True)
A = Symbol('A', real=True)
W = Symbol('W', real=True)
L = Symbol('L', real=True)
M = Symbol('M', real=True)
D1 = Symbol('D1', real=True)
E1 = Symbol('E1', real=True)
A1 = Symbol('A1', real=True)
B_ha = Symbol('B_ha', real=True)
B_hb = Symbol('B_hb', real=True)
Bh = Symbol('Bh', real=True)
B_B = Symbol("B_B", real=True)
λa = Symbol('λa', real=True)
λb = Symbol('λb', real=True)
λ = Symbol('λ', real=True)
α = Symbol('α', real=True)
Ma = Symbol('Ma', real=True)
Mb = Symbol('Mb', real=True)
R_M = Symbol('R_M', real=True)
Ba = Symbol('Ba', real=True)
Bb = Symbol('Bb', real=True)
I = Symbol('I', real=True)
Pf = Symbol('Pf', real=True)
R_I = Symbol('R_I', real=True)
R_L = Symbol('R_L', real=True)
μ_ma = Symbol('μ_ma', real=True)
μ_Ba = Symbol('μ_Ba', real=True)
μ_Bb = Symbol('μ_Bb', real=True)
μ_B = Symbol('μ_B', real=True)
P_B1 = Symbol('P_B1', real=True)
B = Symbol('B', real=True)
P = Symbol('P', real=True)
r = Symbol('r', real=True)

# symbols with default assumptions (complex field)

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





# Lagrangian Multiplier

η = Symbol('η')


# symbols as undefined function

f = symbols("f", cls=Function)
r = symbols('r', cls=Function)




##################### Objective #######################

"""

# Firm obj =

           L
           ⌠
-L⋅R_L + P⋅⎮ r(L) dL
           ⌡
           0
"""


Firm = P * integrate(r(L), (L, 0, L)) - R_L * L

Firm.diff(L).doit()

# F.O.C,   P⋅r(L) - R_L



"""

# MMF obj =

B_H⋅R_B - B_H⋅R_F0

"""

MMF = R_B * B_H - R_F0 * B_H

MMF.diff(B_H)

# F.O.C,   R_B - R_F0


"""

# HH obj



"""

HH = 1 / (2*P) * (2*(R_F0*θ)*F +(R_D*θ)*D + 2*R_E*E + R_E1*E1 \
                  + R_D1*θ*(B1*P_B - E1) + R_D*θ*(D-λ*(D+E))
                  + R_F0*θ*(λ*(D+E)*R_W/P_B) - R_F0*θ*B1)

s_F = W - D - E
s_D1 = B1*P_B -E1

HH = HH.subs(F,s_F).subs(D1,s_D1)

HH.diff(B1)
solve(HH.diff(B1),R_F0, dict=True)

# [{R_F0: P_B⋅R_D1}]

HH.diff(E1)
solve(HH.diff(E1),R_E1, dict=True)

# [{R_E1: R_D1⋅θ}]

HH.diff(E)

"""
                              R_F0⋅R_W⋅θ⋅λ
-R_D⋅θ⋅λ + 2⋅R_E - 2⋅R_F0⋅θ + ────────────
                                  P_B
──────────────────────────────────────────
                   2⋅P
"""

HH.diff(D)

"""
                                    R_F0⋅R_W⋅θ⋅λ
R_D⋅θ⋅(-λ + 1) + R_D⋅θ - 2⋅R_F0⋅θ + ────────────
                                        P_B
────────────────────────────────────────────────
                      2⋅P

"""


"""
Bank

"""


Bank = 1/2 * (R_L*L + R_M**2*M -R_E*E -R_D*D +R_M*A1 -R_D1*D1 \
              -R_E1*E1 +(R_I-R_M)*I - P*integrate(f(I),(I,0,I))) \
       + 1/2 * (R_L*L -R_E*E - R_I*(λ*A*R_W-R_M*M) \
                -R_D*((1-λ)*A-E))


s_L = A - M
s_E = A*α*A / R_E
s_E1 = A1*α*(A+A1)/R_E1
s_D = A - s_E
s_D1 = A1 - s_E1

Bank = Bank.subs(L,s_L).subs(E,s_E).subs(D,s_D).subs(D1,s_D1)

Bank.diff(I).doit()
s_R_I = solve(Bank.diff(I).doit(), R_I)[0]


Bank.diff(M).subs(R_I,s_R_I).simplify()

"""

                                  2
0.5⋅P⋅R_M⋅f(I) - 1.0⋅R_L + 1.0⋅R_M


"""

Bank.diff(A1)
Bank.diff(A)


##################### Equations #######################

# Firm, FOC
feq1 = P - R_L/r(L)

# Money Fund, FOC

meq1 = R_B - R_F0


# HH, FOC

eq1 = R_D1 * θ * P_B - (R_F0 * θ)
eq2 = R_E1 - R_D1 * θ
eq3 = -2*R_F0*θ + 2*R_E - λ*R_D*θ + R_F0 * θ * λ * R_W / P_B
eq4 = -2*R_F0*θ + R_D*θ + (1-λ)*R_D*θ + R_F0 * θ * λ * R_W / P_B

## Bank, FOC

beq1 = R_I - R_M - P*f(I)
beq11 = beq1.subs(f(I),F*I)
beq2 = -2*R_L + R_M*R_M + R_M*R_I
beq3 = R_M - R_D1*(1-(α*(A+A1)+A1*α)/R_E1)-R_E1*(α*(A+A1)+A1*α)/R_E1


s_e0 = 2*α*A
s_e1 = (1+2*λ*R_W)*α*A

beq4 = 2*R_L - 2*R_E*s_e0/R_E - R_D*(1-s_e0/R_E)- R_D1*(A1*α/R_E1) \
       - R_E1*(A1*α/R_E1)

beq5 = A1 - λ*R_W*A




def main():

    ############################# April 10th Experiments ############

    s_1 = solve((eq3, eq4), R_E, R_F0, λ, R_W, P_B)

    """
    ⎡⎧           R_F0⋅R_W⋅λ                  ⎫⎤
    ⎢⎨P_B: ──────────────────────, R_E: R_D⋅θ⎬⎥
    ⎣⎩     R_D⋅λ - 2⋅R_D + 2⋅R_F0            ⎭⎦
    """

    s_1 = solve((eq3, eq4), R_E)

    """
    ⎧     θ⋅(P_B⋅R_D⋅λ + 2⋅P_B⋅R_F0 - R_F0⋅R_W⋅λ)⎫
    ⎨R_E: ───────────────────────────────────────⎬
    ⎩                      2⋅P_B                 ⎭
    """

    s_1 = solve((eq3, eq4), R_E, P_B, dict=True)

    """
    ⎡⎧           R_F0⋅R_W⋅λ                  ⎫⎤
    ⎢⎨P_B: ──────────────────────, R_E: R_D⋅θ⎬⎥
    ⎣⎩     R_D⋅λ - 2⋅R_D + 2⋅R_F0            ⎭⎦
    """

    s_1 = solve((eq3, eq4), R_E, R_F0)

    """
    ⎧                  P_B⋅R_D⋅(-λ + 2)⎫
    ⎨R_E: R_D⋅θ, R_F0: ────────────────⎬
    ⎩                   2⋅P_B - R_W⋅λ  ⎭
    """

    s_1 = solve((eq3, eq4), R_E, R_F0, P_B, dict=True)

    """
    ⎡⎧           R_F0⋅R_W⋅λ                  ⎫⎤
    ⎢⎨P_B: ──────────────────────, R_E: R_D⋅θ⎬⎥
    ⎣⎩     R_D⋅λ - 2⋅R_D + 2⋅R_F0            ⎭⎦
    """

    s_1 = solve((eq3, eq4), R_E, R_F0, P_B, R_D, dict=True)

    """
    ⎡⎧           R_F0⋅R_W⋅θ⋅λ             R_E⎫⎤
    ⎢⎨P_B: ────────────────────────, R_D: ───⎬⎥
    ⎣⎩     R_E⋅λ - 2⋅R_E + 2⋅R_F0⋅θ        θ ⎭⎦
    """

    s_1 = solve((eq3, eq4), R_D, R_E, R_F0, P_B, dict=True)

    """
    ⎡⎧           R_F0⋅R_W⋅θ⋅λ             R_E⎫⎤
    ⎢⎨P_B: ────────────────────────, R_D: ───⎬⎥
    ⎣⎩     R_E⋅λ - 2⋅R_E + 2⋅R_F0⋅θ        θ ⎭⎦
    """

    s_1 = solve((eq3, eq4))

    """
    ⎡⎧           R_F0⋅R_W⋅θ⋅λ             R_E⎫                ⎤
    ⎢⎨P_B: ────────────────────────, R_D: ───⎬, {R_E: 0, θ: 0}⎥
    ⎣⎩     R_E⋅λ - 2⋅R_E + 2⋅R_F0⋅θ        θ ⎭                ⎦
    """

    s_2 = solve((eq1, eq2, eq3), R_E, R_F0, R_E1, P_B)

    """
    ⎡⎧     R_F0       θ⋅(R_D⋅λ - R_D1⋅R_W⋅λ + 2⋅R_F0)              ⎫⎤
    ⎢⎨P_B: ────, R_E: ───────────────────────────────, R_E1: R_D1⋅θ⎬⎥
    ⎣⎩     R_D1                      2                             ⎭⎦
    """


    ############################# April 11th Experiments ############



    # other solvers

    s_1 = solve_poly_system([eq3, eq4], R_E, R_F0)

    """
    ⎡⎛       -P_B⋅R_D⋅(λ - 2) ⎞⎤
    ⎢⎜R_D⋅θ, ─────────────────⎟⎥
    ⎣⎝         2⋅P_B - R_W⋅λ  ⎠⎦
    """

    s_1 = nonlinsolve([eq3, eq4], R_E, P_B)

    """
    ⎧⎛             R_F0⋅R_W⋅λ      ⎞⎫
    ⎨⎜R_D⋅θ, ──────────────────────⎟⎬
    ⎩⎝       R_D⋅λ - 2⋅R_D + 2⋅R_F0⎠⎭
    """

    # undefined function and its derivative

    f = symbols("f", cls=Function)

    f(R_E).diff(R_E)

    """
     d
    ────(f(R_E))
    dR_E
    """


    # taking derivative of defined function

    s_1 = solve((eq3, eq4), dict=True)

    s_1[0][P_B].diff(R_E)

    """
       R_F0⋅R_W⋅θ⋅λ⋅(-λ + 2)
    ───────────────────────────
                              2
    (R_E⋅λ - 2⋅R_E + 2⋅R_F0⋅θ)
    """


    # mimic man's solution

        #Household

    s_R_F0 = solve(eq1,R_F0)[0]
    "P_B⋅R_D1"
    s_R_E1 = solve(eq2, R_E1)[0]
    "R_D1⋅θ"
    s_R_E = solve((eq3, eq4), R_E, P_B, dict=True)[0][R_E]
    "R_D⋅θ"
    s_P_B = simplify(
        solve(eq3.subs(s_R_E, R_E).subs(R_F0, s_R_F0)
              .subs(s_R_E1, R_E1), P_B)[0]
            .subs(R_E, s_R_E).subs(R_E1, s_R_E1))

    """
    -R_D⋅λ + 2⋅R_D + R_D1⋅R_W⋅λ
    ───────────────────────────
               2⋅R_D1
    """
    s2_R_F0 = collect(s_R_F0.subs(P_B,s_P_B), R_D)
    """
        ⎛  λ    ⎞   R_D1⋅R_W⋅λpp = pprint.PrettyPrinter(indent=4)
    R_D⋅⎜- ─ + 1⎟ + ──────────
        ⎝  2    ⎠       2
    """

        # Mutual Fund

    s_R_B = solve(meq1, R_B)[0]
    "R_F0"
    s2_R_B = s_R_B.subs(R_F0,s2_R_F0)
    """
        ⎛  λ    ⎞   R_D1⋅R_W⋅λ
    R_D⋅⎜- ─ + 1⎟ + ──────────
        ⎝  2    ⎠       2

    """

        # Bank


    s_R_I = solve(beq1, R_I)[0]
    "P⋅f(I) + R_M"


    integrate(f(I), (I, 0, I))

    """
    I
    ⌠
    ⎮ f(I) dI
    ⌡
    0
    """
    integrate(f(I), (I, 0, I)).diff(I).doit()
    "f(I)"

    ############################# April 12th Experiments ############

    s_R_L = solve(beq2, R_L)[0].subs(R_I,s_R_I)

    """
    R_M⋅(P⋅f(I) + 2⋅R_M)
    ────────────────────
            2

    """

    s_R_D1 = solve(beq3, R_D1)[0]

    """

    R_E1⋅(A⋅α + 2⋅A₁⋅α - R_M)
    ─────────────────────────
        A⋅α + 2⋅A₁⋅α - R_E1

    """

    s_A1 = solve(beq5, A1)[0]
    s2_R_D1 = s_R_D1.subs(A1,s_A1).subs(s_e1,e1)

    sss = solve((meq1,feq1,eq1,eq2,eq3,eq4,beq1,beq2,beq3,beq4,beq5)
                ,R_D1,R_E1,R_D,R_E,P_B,R_F0,P,R_L,R_I,A,A1
                , dict=True)

    ss2 = solve((meq1, feq1, eq1, eq2, eq3, eq4, beq1, beq2, beq3, beq4, beq5)
                , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I, A, A1, R_B
                , dict=True)

    ss3 = solve((meq1, feq1, eq1, eq2, eq3, eq4, beq1, beq2, beq3, beq4, beq5)
                , dict=True)

    # beq11 = beq1.subs(f(I),F*I)

    ss4 = solve((meq1, feq1, eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4, beq5)
                , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I, R_B , I
                , dict=True)

    ss5 = solve((meq1, feq1, eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4, beq5)
                , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I, R_B
                , dict=True)

    ss5_1 = solve((meq1, feq1, eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I, R_B
                  , dict=True)

    ss5_2 = solve((meq1, feq1, eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I, R_B, I
                  , dict=True)

    ss5_3 = solve(( feq1, eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I, I
                  , dict=True)

    ss5_4 = solve((eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_F0, R_L, R_I, I
                  , dict=True)

    ss5_5 = solve((eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_L, R_I, I
                  , dict=True)

    ss5_6 = solve((eq1, eq2, eq3, eq4, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_L, R_I
                  , dict=True)

    ss5_7 = solve((eq1, eq2, eq3, eq4, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_L
                  , dict=True)

    ss5_8 = solve(( eq2, eq3, eq4, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, R_L
                  , dict=True)

    ss5_9 = solve((eq2, eq3, eq4, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, R_L, R_F0
                  , dict=True)

    ss5_10 = solve((eq2, eq3, eq4, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, R_L
                  , dict=True)

    ss5_11 = solve((eq2, eq3, eq4, beq2, beq3)
                   , R_D1, R_E1, R_D, R_E
                   , dict=True)

    ss5_12 = solve((eq2, eq3, beq2, beq3)
                   , R_D1, R_E1, R_E
                   , dict=True)

    ss5_13 = solve((eq2, eq3, eq4, beq3)
                   , R_D1, R_E1, R_D, R_E
                   , dict=True)

    ss5_14 = solve((eq2, eq3, beq3)
                   , R_D1, R_E1, R_E
                   , dict=True)


    ss5_15 = solve((eq2, eq3, beq2, beq3)
                   , R_D1, R_E1, R_E, R_L
                   , dict=True)

    ss5_16 = solve((eq2, eq3, beq3, beq4)
                   , R_D1, R_E1, R_E, R_D
                   , dict=True)

    ss5_17 = solve((eq2, eq3, beq2, beq3, beq4)
                   , R_D1, R_E1, R_E, R_D, R_L
                   , dict=True)

    ss5_18 = solve((eq2, eq3, eq4, beq2, beq3, beq4)
                   , R_D1, R_E1, R_E, R_D, R_L, R_F0
                   , dict=True)

    ss5_19 = solve((eq2, eq3, eq4, beq3, beq4)
                   , R_D1, R_E1, R_E, R_D, R_F0
                   , dict=True)

    ss5_20 = solve((eq2, eq3, eq4, beq3)
                   , R_D1, R_E1, R_E, R_F0
                   , dict=True)

    ss5_21 = solve(( eq3, eq4,eq1,eq2,beq3,beq4)
                   , R_E, R_F0,P_B,R_E1,R_D1,R_D
                   , dict=True)

    ss5_22 = solve((eq3, eq4, eq1, eq2, beq3, beq2)
                   , R_E, R_F0, P_B, R_E1, R_D1,R_L
                   , dict=True)

    ss5_23 = solve((eq3, eq4, eq1, eq2, beq3, beq2,beq1,feq1)
                   , R_E, R_F0, P_B, R_E1, R_D1, R_L,R_I,P
                   , dict=True)

    ss5_24 = solve((eq3, eq4, eq1, eq2, beq3, beq2, beq1, feq1)
                   , R_E, R_F0, P_B, R_E1, R_D1, R_L, R_I, P
                   , dict=True)




    ss6 = solve((feq1, eq1, eq2, eq3, eq4, beq1, beq2, beq3, beq4, beq5)
                , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I
                , dict=True)

    ss7 = solve((feq1, eq1, eq2, eq3, eq4, beq1, beq2, beq3, beq4)
                , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I
                , dict=True)

    ss8 = solve((meq1, feq1, eq1, eq2, eq3, eq4, beq1, beq2, beq3, beq4, beq5)
                , P, R_L, R_I, R_D1, R_E1, R_D, R_E, P_B, R_F0, A, A1
                , dict=True)


    s_R_D = solve((beq4), R_D)
    s_R_D2 = solve((beq4, eq2), R_D)
    s_R_D3 = solve((beq4, beq2), R_D, R_L)
    s_R_D3 = solve((beq4, beq2, eq2)
                   , R_D, R_L, R_E1
                   , dict=True)
    s_R_D3 = solve((beq4, beq2, eq2)
                   , R_D, R_L, R_E1
                   , dict=True)

    s_R_D4 = solve((beq4, beq2, eq2, eq3)
                   , R_D, R_L, R_E1, R_E
                   , dict=True)

    s_R_D5 = solve((beq4, beq2, eq2, eq3, eq4)
                   , R_D, R_L, R_E1, R_E, R_F0
                   , dict=True)

    s_R_D6 = solve((beq4, beq2, eq2, eq3, eq4, eq1)
                   , R_D, R_L, R_E1, R_E, R_F0, P_B
                   , dict=True)

    s_R_D7 = solve((beq4, beq2, beq1, eq2, eq3, eq4, eq1)
                   , R_D, R_L, R_E1, R_E, R_F0, P_B, R_I
                   , dict=True)

    s_R_D8 = solve((beq4, beq2, eq2, eq3, eq4, eq1, feq1)
                   , R_D, R_L, R_E1, R_E, R_F0, P_B, P
                   , dict=True)

    s_R_D9 = solve((beq4, beq2, beq1, eq2, eq3, eq4, eq1, feq1)
                   , R_D, R_L, R_E1, R_E, R_F0, P_B, R_I, P
                   , dict=True)

    s_R_D10 = solve((beq4, beq3, beq2, beq1, eq2, eq3, eq4, eq1, feq1)
                   , R_D, R_L, R_E1, R_E, R_F0, P_B, R_I, P, R_D1
                   , dict=True)

    s_R_D11 = solve((beq4, beq3)
                    , R_D, R_D1
                    , dict=True)

    s_R_D12 = solve((beq4, beq3, eq2)
                    , R_D, R_D1, R_E1
                    , dict=True)

    s_R_D13 = solve((beq4, beq3, beq2, eq2)
                    , R_D, R_D1, R_E1, R_L
                    , dict=True)

    s_R_D14 = solve((beq4, beq3, eq2, eq3)
                    , R_D, R_D1, R_E1, R_E
                    , dict=True)

    s_R_D15 = solve((beq4, beq3, eq2, eq3, eq4)
                    , R_D, R_D1, R_E1, R_E, R_F0
                    , dict=True)

    s_R_D16 = solve((beq4, beq3, eq2, eq4)
                    , R_D, R_D1, R_E1, R_F0
                    , dict=True)

    s_R_D16 = solve((beq4, beq3, eq2, eq4)
                    , R_D, R_D1, R_E1, R_F0
                    , dict=True)



eqs = [beq4, beq3, beq2, beq1, eq2, eq3, eq4, eq1, feq1]
vars = [R_D, R_L, R_E1, R_E, R_F0, P_B, R_I, P, R_D1]


s_R_D16 = solve((beq4, beq3, beq2, beq1, eq2, eq3, eq1, feq1)
                    , R_D, R_L, R_E1, R_E, P_B, R_I, P, R_D1
                    , dict=True)


s_R_D16 = solve((beq4, beq3, beq2, beq1, eq2, eq3, eq4, eq1)
                    , R_D, R_L, R_E1, R_E, P_B, R_I, R_D1, R_F0
                    , dict=True)

s_R_D17 = solve((beq4, beq3, beq2, beq1, eq2, eq3, feq1)
                    , R_D, R_L, R_E1, R_E, R_I, P, R_D1
                    , dict=True)


s_R_D17 = solve((beq4, beq3, beq2, beq1, eq1, eq3, feq1)
                    , R_D, R_L, R_E, R_I, P, R_D1, P_B
                    , dict=True)


s_R_D17 = solve((beq4, beq3, beq2, beq1, eq1, eq3, eq4)
                    , R_D, R_L, R_E, R_I, R_D1, P_B, R_F0
                    , dict=True)


ss = solve((beq4, beq3, eq1, eq2, eq4)
                    , R_D, R_D1, P_B, R_F0, R_E1
                    , dict=True)

ss = solve((beq4, beq3, eq1, eq2)
                    , R_D, R_E1, R_D1, P_B
                    , dict=True)




def myprint(dict):
    [pprint(v) for k,v in dict.items()]







if __name__ == "__main__":
    main()





"""

********************** RD_1 *********************************

4⋅R_B⋅R_M⋅R_W⋅θ⋅λ⋅f(I) - 4⋅R_B⋅R_M⋅R_W⋅λ⋅f(I) + 2⋅R_B⋅R_M⋅θ⋅f(I) - 2⋅R_B⋅R_M⋅f
──────────────────────────────────────────────────────────────────────────────


                                                                             2
(I) - 8⋅R_B⋅R_W⋅θ⋅λ⋅r(L) + 8⋅R_B⋅R_W⋅λ⋅r(L) - 4⋅R_B⋅θ⋅r(L) + 4⋅R_B⋅r(L) - R_M
──────────────────────────────────────────────────────────────────────────────


        2             2        2             2                      2
⋅R_W⋅θ⋅λ ⋅f(I) - 8⋅R_M ⋅R_W⋅θ⋅λ ⋅r(L) + 2⋅R_M ⋅R_W⋅θ⋅λ⋅f(I) + 16⋅R_M ⋅R_W⋅θ⋅λ⋅
──────────────────────────────────────────────────────────────────────────────


          2      2             2      2             2                    2
r(L) - R_M ⋅R_W⋅λ ⋅f(I) + 8⋅R_M ⋅R_W⋅λ ⋅r(L) + 2⋅R_M ⋅R_W⋅λ⋅f(I) - 16⋅R_M ⋅R_W
──────────────────────────────────────────────────────────────────────────────
                          ⎛     2    2        2  2          2
      (R_M⋅f(I) - 2⋅r(L))⋅⎝2⋅R_W ⋅θ⋅λ  - 2⋅R_W ⋅λ  - R_W⋅θ⋅λ  + 3⋅R_W⋅θ⋅λ - R_
               2                 2                 2               2
⋅λ⋅r(L) - 4⋅R_M ⋅θ⋅λ⋅f(I) - 4⋅R_M ⋅θ⋅λ⋅r(L) + 8⋅R_M ⋅θ⋅f(I) + 8⋅R_M ⋅θ⋅r(L) +
──────────────────────────────────────────────────────────────────────────────
   2                                ⎞
W⋅λ  + R_W⋅λ - 4⋅θ⋅λ + 8⋅θ + 2⋅λ - 4⎠
     2               2               2             2                     2
2⋅R_M ⋅λ⋅f(I) + 4⋅R_M ⋅λ⋅r(L) - 4⋅R_M ⋅f(I) - 8⋅R_M ⋅r(L) + 2⋅R_M⋅R_W⋅θ⋅λ ⋅r(L
──────────────────────────────────────────────────────────────────────────────


                                    2
) - 4⋅R_M⋅R_W⋅θ⋅λ⋅r(L) + 2⋅R_M⋅R_W⋅λ ⋅r(L) - 4⋅R_M⋅R_W⋅λ⋅r(L) + 8⋅R_M⋅θ⋅λ⋅r(L)
──────────────────────────────────────────────────────────────────────────────



 - 16⋅R_M⋅θ⋅r(L) - 4⋅R_M⋅λ⋅r(L) + 8⋅R_M⋅r(L)
────────────────────────────────────────────



********************** RE_1 *********************************

  ⎛
θ⋅⎝4⋅R_B⋅R_M⋅R_W⋅θ⋅λ⋅f(I) - 4⋅R_B⋅R_M⋅R_W⋅λ⋅f(I) + 2⋅R_B⋅R_M⋅θ⋅f(I) - 2⋅R_B⋅R_
──────────────────────────────────────────────────────────────────────────────



M⋅f(I) - 8⋅R_B⋅R_W⋅θ⋅λ⋅r(L) + 8⋅R_B⋅R_W⋅λ⋅r(L) - 4⋅R_B⋅θ⋅r(L) + 4⋅R_B⋅r(L) - R
──────────────────────────────────────────────────────────────────────────────


  2        2             2        2             2                      2
_M ⋅R_W⋅θ⋅λ ⋅f(I) - 8⋅R_M ⋅R_W⋅θ⋅λ ⋅r(L) + 2⋅R_M ⋅R_W⋅θ⋅λ⋅f(I) + 16⋅R_M ⋅R_W⋅θ
──────────────────────────────────────────────────────────────────────────────


             2      2             2      2             2                    2
⋅λ⋅r(L) - R_M ⋅R_W⋅λ ⋅f(I) + 8⋅R_M ⋅R_W⋅λ ⋅r(L) + 2⋅R_M ⋅R_W⋅λ⋅f(I) - 16⋅R_M ⋅
──────────────────────────────────────────────────────────────────────────────
                            ⎛     2    2        2  2          2
        (R_M⋅f(I) - 2⋅r(L))⋅⎝2⋅R_W ⋅θ⋅λ  - 2⋅R_W ⋅λ  - R_W⋅θ⋅λ  + 3⋅R_W⋅θ⋅λ -
                  2                 2                 2               2
R_W⋅λ⋅r(L) - 4⋅R_M ⋅θ⋅λ⋅f(I) - 4⋅R_M ⋅θ⋅λ⋅r(L) + 8⋅R_M ⋅θ⋅f(I) + 8⋅R_M ⋅θ⋅r(L)
──────────────────────────────────────────────────────────────────────────────
     2                                ⎞
R_W⋅λ  + R_W⋅λ - 4⋅θ⋅λ + 8⋅θ + 2⋅λ - 4⎠
        2               2               2             2                     2
 + 2⋅R_M ⋅λ⋅f(I) + 4⋅R_M ⋅λ⋅r(L) - 4⋅R_M ⋅f(I) - 8⋅R_M ⋅r(L) + 2⋅R_M⋅R_W⋅θ⋅λ ⋅
──────────────────────────────────────────────────────────────────────────────


                                       2
r(L) - 4⋅R_M⋅R_W⋅θ⋅λ⋅r(L) + 2⋅R_M⋅R_W⋅λ ⋅r(L) - 4⋅R_M⋅R_W⋅λ⋅r(L) + 8⋅R_M⋅θ⋅λ⋅r
──────────────────────────────────────────────────────────────────────────────


                                               ⎞
(L) - 16⋅R_M⋅θ⋅r(L) - 4⋅R_M⋅λ⋅r(L) + 8⋅R_M⋅r(L)⎠
────────────────────────────────────────────────


********************** RD   *********************************

2⋅R_B⋅R_M⋅R_W⋅θ⋅λ⋅f(I) + 2⋅R_B⋅R_M⋅R_W⋅λ⋅f(I) + 8⋅R_B⋅R_M⋅θ⋅f(I) - 4⋅R_B⋅R_M⋅f
──────────────────────────────────────────────────────────────────────────────



(I) - 4⋅R_B⋅R_W⋅θ⋅λ⋅r(L) - 4⋅R_B⋅R_W⋅λ⋅r(L) - 16⋅R_B⋅θ⋅r(L) + 8⋅R_B⋅r(L) - R_M
──────────────────────────────────────────────────────────────────────────────

                                                                      (R_M⋅f(I
2    2    2             2    2    2           2    2  2             2    2  2
 ⋅R_W ⋅θ⋅λ ⋅f(I) - 8⋅R_M ⋅R_W ⋅θ⋅λ ⋅r(L) - R_M ⋅R_W ⋅λ ⋅f(I) + 8⋅R_M ⋅R_W ⋅λ ⋅
──────────────────────────────────────────────────────────────────────────────
            ⎛     2    2        2  2          2                    2
) - 2⋅r(L))⋅⎝2⋅R_W ⋅θ⋅λ  - 2⋅R_W ⋅λ  - R_W⋅θ⋅λ  + 3⋅R_W⋅θ⋅λ - R_W⋅λ  + R_W⋅λ -
            2                     2                     2                   2
r(L) - 4⋅R_M ⋅R_W⋅θ⋅λ⋅f(I) - 4⋅R_M ⋅R_W⋅θ⋅λ⋅r(L) + 2⋅R_M ⋅R_W⋅λ⋅f(I) + 4⋅R_M ⋅
──────────────────────────────────────────────────────────────────────────────
                      ⎞
 4⋅θ⋅λ + 8⋅θ + 2⋅λ - 4⎠
                      2    2                 2  2
R_W⋅λ⋅r(L) + 2⋅R_M⋅R_W ⋅θ⋅λ ⋅r(L) + 2⋅R_M⋅R_W ⋅λ ⋅r(L) + 8⋅R_M⋅R_W⋅θ⋅λ⋅r(L) -
──────────────────────────────────────────────────────────────────────────────



4⋅R_M⋅R_W⋅λ⋅r(L)
────────────────


********************** RE   *********************************

  ⎛
θ⋅⎝2⋅R_B⋅R_M⋅R_W⋅θ⋅λ⋅f(I) + 2⋅R_B⋅R_M⋅R_W⋅λ⋅f(I) + 8⋅R_B⋅R_M⋅θ⋅f(I) - 4⋅R_B⋅R_
──────────────────────────────────────────────────────────────────────────────



M⋅f(I) - 4⋅R_B⋅R_W⋅θ⋅λ⋅r(L) - 4⋅R_B⋅R_W⋅λ⋅r(L) - 16⋅R_B⋅θ⋅r(L) + 8⋅R_B⋅r(L) -
──────────────────────────────────────────────────────────────────────────────

                                                                        (R_M⋅f
   2    2    2             2    2    2           2    2  2             2    2
R_M ⋅R_W ⋅θ⋅λ ⋅f(I) - 8⋅R_M ⋅R_W ⋅θ⋅λ ⋅r(L) - R_M ⋅R_W ⋅λ ⋅f(I) + 8⋅R_M ⋅R_W ⋅
──────────────────────────────────────────────────────────────────────────────
              ⎛     2    2        2  2          2                    2
(I) - 2⋅r(L))⋅⎝2⋅R_W ⋅θ⋅λ  - 2⋅R_W ⋅λ  - R_W⋅θ⋅λ  + 3⋅R_W⋅θ⋅λ - R_W⋅λ  + R_W⋅λ
 2             2                     2                     2
λ ⋅r(L) - 4⋅R_M ⋅R_W⋅θ⋅λ⋅f(I) - 4⋅R_M ⋅R_W⋅θ⋅λ⋅r(L) + 2⋅R_M ⋅R_W⋅λ⋅f(I) + 4⋅R_
──────────────────────────────────────────────────────────────────────────────
                        ⎞
 - 4⋅θ⋅λ + 8⋅θ + 2⋅λ - 4⎠
 2                       2    2                 2  2
M ⋅R_W⋅λ⋅r(L) + 2⋅R_M⋅R_W ⋅θ⋅λ ⋅r(L) + 2⋅R_M⋅R_W ⋅λ ⋅r(L) + 8⋅R_M⋅R_W⋅θ⋅λ⋅r(L)
──────────────────────────────────────────────────────────────────────────────


                   ⎞
 - 4⋅R_M⋅R_W⋅λ⋅r(L)⎠
────────────────────


********************** P_B *********************************



──────────────────────────────────────────────────────────────────────────────

4⋅R_B⋅R_M⋅R_W⋅θ⋅λ⋅f(I) - 4⋅R_B⋅R_M⋅R_W⋅λ⋅f(I) + 2⋅R_B⋅R_M⋅θ⋅f(I) - 2⋅R_B⋅R_M⋅f


──────────────────────────────────────────────────────────────────────────────
                                                                             2
(I) - 8⋅R_B⋅R_W⋅θ⋅λ⋅r(L) + 8⋅R_B⋅R_W⋅λ⋅r(L) - 4⋅R_B⋅θ⋅r(L) + 4⋅R_B⋅r(L) - R_M


──────────────────────────────────────────────────────────────────────────────
        2             2        2             2                      2
⋅R_W⋅θ⋅λ ⋅f(I) - 8⋅R_M ⋅R_W⋅θ⋅λ ⋅r(L) + 2⋅R_M ⋅R_W⋅θ⋅λ⋅f(I) + 16⋅R_M ⋅R_W⋅θ⋅λ⋅
                            ⎛     2    2        2  2          2
    R_B⋅(R_M⋅f(I) - 2⋅r(L))⋅⎝2⋅R_W ⋅θ⋅λ  - 2⋅R_W ⋅λ  - R_W⋅θ⋅λ  + 3⋅R_W⋅θ⋅λ -
──────────────────────────────────────────────────────────────────────────────
          2      2             2      2             2                    2
r(L) - R_M ⋅R_W⋅λ ⋅f(I) + 8⋅R_M ⋅R_W⋅λ ⋅r(L) + 2⋅R_M ⋅R_W⋅λ⋅f(I) - 16⋅R_M ⋅R_W
     2                                ⎞
R_W⋅λ  + R_W⋅λ - 4⋅θ⋅λ + 8⋅θ + 2⋅λ - 4⎠
──────────────────────────────────────────────────────────────────────────────
               2                 2                 2               2
⋅λ⋅r(L) - 4⋅R_M ⋅θ⋅λ⋅f(I) - 4⋅R_M ⋅θ⋅λ⋅r(L) + 8⋅R_M ⋅θ⋅f(I) + 8⋅R_M ⋅θ⋅r(L) +


──────────────────────────────────────────────────────────────────────────────
     2               2               2             2                     2
2⋅R_M ⋅λ⋅f(I) + 4⋅R_M ⋅λ⋅r(L) - 4⋅R_M ⋅f(I) - 8⋅R_M ⋅r(L) + 2⋅R_M⋅R_W⋅θ⋅λ ⋅r(L


──────────────────────────────────────────────────────────────────────────────
                                    2
) - 4⋅R_M⋅R_W⋅θ⋅λ⋅r(L) + 2⋅R_M⋅R_W⋅λ ⋅r(L) - 4⋅R_M⋅R_W⋅λ⋅r(L) + 8⋅R_M⋅θ⋅λ⋅r(L)


────────────────────────────────────────────

 - 16⋅R_M⋅θ⋅r(L) - 4⋅R_M⋅λ⋅r(L) + 8⋅R_M⋅r(L)


********************** R_F0 *********************************

R_B


********************** P    *********************************

           2
     -2⋅R_M
─────────────────
R_M⋅f(I) - 2⋅r(L)


********************** R_L *********************************

        2
  -2⋅R_M ⋅r(L)
─────────────────
R_M⋅f(I) - 2⋅r(L)


********************** R_I *********************************

-R_M⋅(R_M⋅f(I) + 2⋅r(L))
─────────────────────────
    R_M⋅f(I) - 2⋅r(L)


********************** A  *********************************

         ⎛                                 2                   2
      -θ⋅⎝2⋅R_B⋅R_M⋅f(I) - 4⋅R_B⋅r(L) - R_M ⋅R_W⋅λ⋅f(I) - 4⋅R_M ⋅λ⋅r(L) + 8⋅R_
──────────────────────────────────────────────────────────────────────────────
                      ⎛     2    2        2  2          2                    2
α⋅(R_M⋅f(I) - 2⋅r(L))⋅⎝2⋅R_W ⋅θ⋅λ  - 2⋅R_W ⋅λ  - R_W⋅θ⋅λ  + 3⋅R_W⋅θ⋅λ - R_W⋅λ
 2                        ⎞
M ⋅r(L) + 2⋅R_M⋅R_W⋅λ⋅r(L)⎠
─────────────────────────────────
                                ⎞
 + R_W⋅λ - 4⋅θ⋅λ + 8⋅θ + 2⋅λ - 4⎠


********************** A1   *********************************



──────────────────────────────────────────────────────────────────────────────
  ⎛         2    2                 2  2                   2
α⋅⎝2⋅R_M⋅R_W ⋅θ⋅λ ⋅f(I) - 2⋅R_M⋅R_W ⋅λ ⋅f(I) - R_M⋅R_W⋅θ⋅λ ⋅f(I) + 3⋅R_M⋅R_W⋅θ
                                          ⎛                                 2
                                 -R_W⋅θ⋅λ⋅⎝2⋅R_B⋅R_M⋅f(I) - 4⋅R_B⋅r(L) - R_M ⋅
──────────────────────────────────────────────────────────────────────────────
                   2
⋅λ⋅f(I) - R_M⋅R_W⋅λ ⋅f(I) + R_M⋅R_W⋅λ⋅f(I) - 4⋅R_M⋅θ⋅λ⋅f(I) + 8⋅R_M⋅θ⋅f(I) + 2
                  2               2                        ⎞
R_W⋅λ⋅f(I) - 4⋅R_M ⋅λ⋅r(L) + 8⋅R_M ⋅r(L) + 2⋅R_M⋅R_W⋅λ⋅r(L)⎠
──────────────────────────────────────────────────────────────────────────────
                                2    2             2  2                 2
⋅R_M⋅λ⋅f(I) - 4⋅R_M⋅f(I) - 4⋅R_W ⋅θ⋅λ ⋅r(L) + 4⋅R_W ⋅λ ⋅r(L) + 2⋅R_W⋅θ⋅λ ⋅r(L)


──────────────────────────────────────────────────────────────────────────────
                           2
 - 6⋅R_W⋅θ⋅λ⋅r(L) + 2⋅R_W⋅λ ⋅r(L) - 2⋅R_W⋅λ⋅r(L) + 8⋅θ⋅λ⋅r(L) - 16⋅θ⋅r(L) - 4⋅


────────────────
               ⎞
λ⋅r(L) + 8⋅r(L)⎠


********************** RD_1 *********************************

********************** RD_1 *********************************

********************** RD_1 *********************************

********************** RD_1 *********************************

********************** RD_1 *********************************

********************** RD_1 *********************************

********************** RD_1 *********************************

********************** RD_1 *********************************

********************** RD_1 *********************************




"""



"""

⎡⎧
⎢⎪     ⎛
⎢⎨   2⋅⎝A⋅R_W⋅α⋅θ⋅λ - A⋅R_W⋅α⋅λ - 4⋅A⋅α⋅θ⋅λ + 8⋅A⋅α⋅θ + 2⋅A⋅α⋅λ - 4⋅A⋅α + 2⋅A₁
⎢⎪I: ─────────────────────────────────────────────────────────────────────────
⎣⎩                F⋅R_M⋅(A⋅R_W⋅α⋅θ⋅λ - A⋅R_W⋅α⋅λ - 4⋅A⋅α⋅θ⋅λ + 8⋅A⋅α⋅θ + 2⋅A⋅α


⋅R_W⋅α⋅θ⋅λ - 2⋅A₁⋅R_W⋅α⋅λ - A₁⋅α⋅θ⋅λ + 2⋅A₁⋅α⋅θ - A₁⋅α⋅λ + 2⋅A₁⋅α + 2⋅R_F0⋅θ +
──────────────────────────────────────────────────────────────────────────────
⋅λ - 4⋅A⋅α + 2⋅A₁⋅R_W⋅α⋅θ⋅λ - 2⋅A₁⋅R_W⋅α⋅λ - A₁⋅α⋅θ⋅λ + 2⋅A₁⋅α⋅θ - A₁⋅α⋅λ + 2⋅
                                                 A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
      2            2                ⎞          - ─────────── + ───────── + 2⋅A
 2⋅R_M ⋅θ⋅λ - 4⋅R_M ⋅θ - R_M⋅R_W⋅θ⋅λ⎠⋅r(L)            2            2
──────────────────────────────────────────, P: ───────────────────────────────
A₁⋅α + 2⋅R_F0⋅θ - R_M⋅R_W⋅θ⋅λ)
                                                               A₁⋅α⋅θ⋅λ
⋅α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅A⋅α - A₁⋅R_W⋅α⋅θ⋅λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α
                                                                  2
──────────────────────────────────────────────────────────────────────────────
                                     θ⋅(λ - 2)⋅r(L)
     A₁⋅α⋅λ                   R_M⋅R_W⋅θ⋅λ
⋅θ + ────── - A₁⋅α - R_F0⋅θ + ───────────
       2                           2                             R_F0⋅θ
─────────────────────────────────────────, P_B: ──────────────────────────────
                                                -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁


                            -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R
──────────, R_B: R_F0, R_D: ──────────────────────────────────────────────────
⋅α + R_M⋅θ                                                      θ⋅(λ - 2)


_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ        A⋅α + 2⋅A₁⋅α - θ⋅(A⋅α + 2⋅A₁⋅α - R_M)
───────────────────────────────, R_D1: ─────────────────────────────────────,
                                                         θ


     -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M
R_E: ─────────────────────────────────────────────────────────────────────────
                                           λ - 2


⋅R_W⋅θ⋅λ                                                       -A⋅R_W⋅α⋅θ⋅λ +
────────, R_E1: -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ, R_I: ───────────────



A⋅R_W⋅α⋅λ + 4⋅A⋅α⋅θ⋅λ - 8⋅A⋅α⋅θ - 2⋅A⋅α⋅λ + 4⋅A⋅α - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅
──────────────────────────────────────────────────────────────────────────────
                                                                      R_M⋅θ⋅(λ

                                                            2            2
α⋅λ + A₁⋅α⋅θ⋅λ - 2⋅A₁⋅α⋅θ + A₁⋅α⋅λ - 2⋅A₁⋅α - 2⋅R_F0⋅θ - R_M ⋅θ⋅λ + 2⋅R_M ⋅θ +
──────────────────────────────────────────────────────────────────────────────
 - 2)
                     A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
                   - ─────────── + ───────── + 2⋅A⋅α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2
 R_M⋅R_W⋅θ⋅λ              2            2
────────────, R_L: ───────────────────────────────────────────────────────────

                                   A₁⋅α⋅θ⋅λ            A₁⋅α⋅λ
⋅A⋅α - A₁⋅R_W⋅α⋅θ⋅λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α⋅θ + ────── - A₁⋅α - R_F0⋅θ
                                      2                  2
──────────────────────────────────────────────────────────────────────────────
           θ⋅(λ - 2)
  R_M⋅R_W⋅θ⋅λ⎫⎤
+ ───────────⎪⎥
       2     ⎬⎥
─────────────⎪⎥
             ⎭⎦



"""

"""

⎧
⎪     ⎛
⎨   2⋅⎝A⋅R_W⋅α⋅θ⋅λ - A⋅R_W⋅α⋅λ - 4⋅A⋅α⋅θ⋅λ + 8⋅A⋅α⋅θ + 2⋅A⋅α⋅λ - 4⋅A⋅α + 2⋅A₁⋅
⎪I: ──────────────────────────────────────────────────────────────────────────
⎩                F⋅R_M⋅(A⋅R_W⋅α⋅θ⋅λ - A⋅R_W⋅α⋅λ - 4⋅A⋅α⋅θ⋅λ + 8⋅A⋅α⋅θ + 2⋅A⋅α⋅


R_W⋅α⋅θ⋅λ - 2⋅A₁⋅R_W⋅α⋅λ - A₁⋅α⋅θ⋅λ + 2⋅A₁⋅α⋅θ - A₁⋅α⋅λ + 2⋅A₁⋅α + 2⋅R_F0⋅θ +
──────────────────────────────────────────────────────────────────────────────
λ - 4⋅A⋅α + 2⋅A₁⋅R_W⋅α⋅θ⋅λ - 2⋅A₁⋅R_W⋅α⋅λ - A₁⋅α⋅θ⋅λ + 2⋅A₁⋅α⋅θ - A₁⋅α⋅λ + 2⋅A
                                                A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
     2            2                ⎞          - ─────────── + ───────── + 2⋅A⋅
2⋅R_M ⋅θ⋅λ - 4⋅R_M ⋅θ - R_M⋅R_W⋅θ⋅λ⎠⋅r(L)            2            2
─────────────────────────────────────────, P: ────────────────────────────────
₁⋅α + 2⋅R_F0⋅θ - R_M⋅R_W⋅θ⋅λ)
                                                              A₁⋅α⋅θ⋅λ
α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅A⋅α - A₁⋅R_W⋅α⋅θ⋅λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α⋅
                                                                 2
──────────────────────────────────────────────────────────────────────────────
                                    θ⋅(λ - 2)⋅r(L)
    A₁⋅α⋅λ                   R_M⋅R_W⋅θ⋅λ
θ + ────── - A₁⋅α - R_F0⋅θ + ───────────
      2                           2                             R_F0⋅θ
────────────────────────────────────────, P_B: ───────────────────────────────
                                               -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅


                           -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_
─────────, R_B: R_F0, R_D: ───────────────────────────────────────────────────
α + R_M⋅θ                                                      θ⋅(λ - 2)


W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ        A⋅α + 2⋅A₁⋅α - θ⋅(A⋅α + 2⋅A₁⋅α - R_M)
──────────────────────────────, R_D1: ─────────────────────────────────────, R
                                                        θ


    -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅
_E: ──────────────────────────────────────────────────────────────────────────
                                          λ - 2


R_W⋅θ⋅λ                                                       -A⋅R_W⋅α⋅θ⋅λ + A
───────, R_E1: -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ, R_I: ────────────────



⋅R_W⋅α⋅λ + 4⋅A⋅α⋅θ⋅λ - 8⋅A⋅α⋅θ - 2⋅A⋅α⋅λ + 4⋅A⋅α - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α
──────────────────────────────────────────────────────────────────────────────
                                                                     R_M⋅θ⋅(λ

                                                           2            2
⋅λ + A₁⋅α⋅θ⋅λ - 2⋅A₁⋅α⋅θ + A₁⋅α⋅λ - 2⋅A₁⋅α - 2⋅R_F0⋅θ - R_M ⋅θ⋅λ + 2⋅R_M ⋅θ +
──────────────────────────────────────────────────────────────────────────────
- 2)
                    A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
                  - ─────────── + ───────── + 2⋅A⋅α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅
R_M⋅R_W⋅θ⋅λ              2            2
───────────, R_L: ────────────────────────────────────────────────────────────

                                  A₁⋅α⋅θ⋅λ            A₁⋅α⋅λ
A⋅α - A₁⋅R_W⋅α⋅θ⋅λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α⋅θ + ────── - A₁⋅α - R_F0⋅θ +
                                     2                  2
──────────────────────────────────────────────────────────────────────────────
          θ⋅(λ - 2)
 R_M⋅R_W⋅θ⋅λ⎫
 ───────────⎪
      2     ⎬
────────────⎪
            ⎭


"""



"""


;;;;;;;;;;;;;;;;;;; ss5_3


ss5_3 = solve(( feq1, eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_F0, P, R_L, R_I, I
                  , dict=True)



⎡⎧
⎢⎪     ⎛
⎢⎨   2⋅⎝A⋅R_W⋅α⋅θ⋅λ - A⋅R_W⋅α⋅λ - 4⋅A⋅α⋅θ⋅λ + 8⋅A⋅α⋅θ + 2⋅A⋅α⋅λ - 4⋅A⋅α + 2⋅A₁
⎢⎪I: ─────────────────────────────────────────────────────────────────────────
⎣⎩                F⋅R_M⋅(A⋅R_W⋅α⋅θ⋅λ - A⋅R_W⋅α⋅λ - 4⋅A⋅α⋅θ⋅λ + 8⋅A⋅α⋅θ + 2⋅A⋅α


⋅R_W⋅α⋅θ⋅λ - 2⋅A₁⋅R_W⋅α⋅λ - A₁⋅α⋅θ⋅λ + 2⋅A₁⋅α⋅θ - A₁⋅α⋅λ + 2⋅A₁⋅α + 2⋅R_F0⋅θ +
──────────────────────────────────────────────────────────────────────────────
⋅λ - 4⋅A⋅α + 2⋅A₁⋅R_W⋅α⋅θ⋅λ - 2⋅A₁⋅R_W⋅α⋅λ - A₁⋅α⋅θ⋅λ + 2⋅A₁⋅α⋅θ - A₁⋅α⋅λ + 2⋅
                                                 A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
      2            2                ⎞          - ─────────── + ───────── + 2⋅A
 2⋅R_M ⋅θ⋅λ - 4⋅R_M ⋅θ - R_M⋅R_W⋅θ⋅λ⎠⋅r(L)            2            2
──────────────────────────────────────────, P: ───────────────────────────────
A₁⋅α + 2⋅R_F0⋅θ - R_M⋅R_W⋅θ⋅λ)
                                                               A₁⋅α⋅θ⋅λ
⋅α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅A⋅α - A₁⋅R_W⋅α⋅θ⋅λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α
                                                                  2
──────────────────────────────────────────────────────────────────────────────
                                     θ⋅(λ - 2)⋅r(L)
     A₁⋅α⋅λ                   R_M⋅R_W⋅θ⋅λ
⋅θ + ────── - A₁⋅α - R_F0⋅θ + ───────────
       2                           2                             R_F0⋅θ
─────────────────────────────────────────, P_B: ──────────────────────────────
                                                -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁


                 -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅
──────────, R_D: ─────────────────────────────────────────────────────────────
⋅α + R_M⋅θ                                           θ⋅(λ - 2)


R_F0⋅θ + R_M⋅R_W⋅θ⋅λ        A⋅α + 2⋅A₁⋅α - θ⋅(A⋅α + 2⋅A₁⋅α - R_M)       -A⋅R_W
────────────────────, R_D1: ─────────────────────────────────────, R_E: ──────
                                              θ


⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ
───────────────────────────────────────────────────────────────────────────, R
                                λ - 2


                                                    -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ +
_E1: -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ, R_I: ──────────────────────────



 4⋅A⋅α⋅θ⋅λ - 8⋅A⋅α⋅θ - 2⋅A⋅α⋅λ + 4⋅A⋅α - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ + A₁⋅α⋅
──────────────────────────────────────────────────────────────────────────────
                                                           R_M⋅θ⋅(λ - 2)

                                                 2            2
θ⋅λ - 2⋅A₁⋅α⋅θ + A₁⋅α⋅λ - 2⋅A₁⋅α - 2⋅R_F0⋅θ - R_M ⋅θ⋅λ + 2⋅R_M ⋅θ + R_M⋅R_W⋅θ⋅
──────────────────────────────────────────────────────────────────────────────

          A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
        - ─────────── + ───────── + 2⋅A⋅α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅A⋅α - A₁⋅R
λ              2            2
─, R_L: ──────────────────────────────────────────────────────────────────────

                        A₁⋅α⋅θ⋅λ            A₁⋅α⋅λ                   R_M⋅R_W⋅θ
_W⋅α⋅θ⋅λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α⋅θ + ────── - A₁⋅α - R_F0⋅θ + ─────────
                           2                  2                           2
──────────────────────────────────────────────────────────────────────────────
θ⋅(λ - 2)
⋅λ⎫⎤
──⎪⎥
  ⎬⎥
──⎪⎥
  ⎭⎦





"""



"""


;;;;;;;;;;;;;;;;;;;;;;;;; ss5_4


    ss5_4 = solve((eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_F0, R_L, R_I, I
                  , dict=True)


⎡⎧
⎢⎪
⎢⎨   -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ + 4⋅A⋅α⋅θ⋅λ - 8⋅A⋅α⋅θ - 2⋅A⋅α⋅λ + 4⋅A⋅α - 2⋅A₁⋅R
⎢⎪I: ─────────────────────────────────────────────────────────────────────────
⎣⎩


_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ + A₁⋅α⋅θ⋅λ - 2⋅A₁⋅α⋅θ + A₁⋅α⋅λ - 2⋅A₁⋅α - 2⋅R_F0⋅θ - 2
──────────────────────────────────────────────────────────────────────────────
           F⋅P⋅R_M⋅θ⋅(λ - 2)

    2            2
⋅R_M ⋅θ⋅λ + 4⋅R_M ⋅θ + R_M⋅R_W⋅θ⋅λ                        R_F0⋅θ
──────────────────────────────────, P_B: ─────────────────────────────────────
                                         -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_


          -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ
───, R_D: ────────────────────────────────────────────────────────────────────
M⋅θ                                           θ⋅(λ - 2)


+ R_M⋅R_W⋅θ⋅λ        A⋅α + 2⋅A₁⋅α - θ⋅(A⋅α + 2⋅A₁⋅α - R_M)       -A⋅R_W⋅α⋅θ⋅λ
─────────────, R_D1: ─────────────────────────────────────, R_E: ─────────────
                                       θ


+ A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ
────────────────────────────────────────────────────────────────────, R_E1: -A
                         λ - 2


                                             -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ + 4⋅A⋅α⋅
⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ, R_I: ─────────────────────────────────



θ⋅λ - 8⋅A⋅α⋅θ - 2⋅A⋅α⋅λ + 4⋅A⋅α - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ + A₁⋅α⋅θ⋅λ - 2
──────────────────────────────────────────────────────────────────────────────
                                                    R_M⋅θ⋅(λ - 2)

                                          2            2
⋅A₁⋅α⋅θ + A₁⋅α⋅λ - 2⋅A₁⋅α - 2⋅R_F0⋅θ - R_M ⋅θ⋅λ + 2⋅R_M ⋅θ + R_M⋅R_W⋅θ⋅λ
────────────────────────────────────────────────────────────────────────, R_L:

   A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
 - ─────────── + ───────── + 2⋅A⋅α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅A⋅α - A₁⋅R_W⋅α⋅θ⋅
        2            2
 ─────────────────────────────────────────────────────────────────────────────
                                                                       θ⋅(λ -
                 A₁⋅α⋅θ⋅λ            A₁⋅α⋅λ                   R_M⋅R_W⋅θ⋅λ⎫⎤
λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α⋅θ + ────── - A₁⋅α - R_F0⋅θ + ───────────⎪⎥
                    2                  2                           2     ⎬⎥
─────────────────────────────────────────────────────────────────────────⎪⎥
2)




"""



"""

;;;;;;;;;;;;;;;;;;;;;;;;; ss5_5


ss5_5 = solve((eq1, eq2, eq3, eq4, beq11, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_L, R_I, I
                  , dict=True)

⎡⎧
⎢⎪
⎢⎨   -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ + 4⋅A⋅α⋅θ⋅λ - 8⋅A⋅α⋅θ - 2⋅A⋅α⋅λ + 4⋅A⋅α - 2⋅A₁⋅R
⎢⎪I: ─────────────────────────────────────────────────────────────────────────
⎣⎩


_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ + A₁⋅α⋅θ⋅λ - 2⋅A₁⋅α⋅θ + A₁⋅α⋅λ - 2⋅A₁⋅α - 2⋅R_F0⋅θ - 2
──────────────────────────────────────────────────────────────────────────────
           F⋅P⋅R_M⋅θ⋅(λ - 2)

    2            2
⋅R_M ⋅θ⋅λ + 4⋅R_M ⋅θ + R_M⋅R_W⋅θ⋅λ                        R_F0⋅θ
──────────────────────────────────, P_B: ─────────────────────────────────────
                                         -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_


          -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ
───, R_D: ────────────────────────────────────────────────────────────────────
M⋅θ                                           θ⋅(λ - 2)


+ R_M⋅R_W⋅θ⋅λ        A⋅α + 2⋅A₁⋅α - θ⋅(A⋅α + 2⋅A₁⋅α - R_M)       -A⋅R_W⋅α⋅θ⋅λ
─────────────, R_D1: ─────────────────────────────────────, R_E: ─────────────
                                       θ


+ A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ
────────────────────────────────────────────────────────────────────, R_E1: -A
                         λ - 2


                                             -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ + 4⋅A⋅α⋅
⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ, R_I: ─────────────────────────────────



θ⋅λ - 8⋅A⋅α⋅θ - 2⋅A⋅α⋅λ + 4⋅A⋅α - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ + A₁⋅α⋅θ⋅λ - 2
──────────────────────────────────────────────────────────────────────────────
                                                    R_M⋅θ⋅(λ - 2)

                                          2            2
⋅A₁⋅α⋅θ + A₁⋅α⋅λ - 2⋅A₁⋅α - 2⋅R_F0⋅θ - R_M ⋅θ⋅λ + 2⋅R_M ⋅θ + R_M⋅R_W⋅θ⋅λ
────────────────────────────────────────────────────────────────────────, R_L:

   A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
 - ─────────── + ───────── + 2⋅A⋅α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅A⋅α - A₁⋅R_W⋅α⋅θ⋅
        2            2
 ─────────────────────────────────────────────────────────────────────────────
                                                                       θ⋅(λ -
                 A₁⋅α⋅θ⋅λ            A₁⋅α⋅λ                   R_M⋅R_W⋅θ⋅λ⎫⎤
λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α⋅θ + ────── - A₁⋅α - R_F0⋅θ + ───────────⎪⎥
                    2                  2                           2     ⎬⎥
─────────────────────────────────────────────────────────────────────────⎪⎥
2)


"""



"""

;;;;;;;;;;;;;;;;;;;;;;;;;;; ss5_6


ss5_6 = solve((eq1, eq2, eq3, eq4, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_L, R_I
                  , dict=True)



⎡⎧
⎢⎪
⎢⎨                      R_F0⋅θ                        -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ
⎢⎪P_B: ────────────────────────────────────────, R_D: ────────────────────────
⎣⎩     -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ


 - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ        A⋅α + 2⋅A₁⋅α
─────────────────────────────────────────────────────────, R_D1: ─────────────
            θ⋅(λ - 2)


- θ⋅(A⋅α + 2⋅A₁⋅α - R_M)       -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A
────────────────────────, R_E: ───────────────────────────────────────────────
     θ                                                               λ - 2


₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ
──────────────────────────────────, R_E1: -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R



           -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ + 4⋅A⋅α⋅θ⋅λ - 8⋅A⋅α⋅θ - 2⋅A⋅α⋅λ + 4⋅A⋅α -
_M⋅θ, R_I: ───────────────────────────────────────────────────────────────────



2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ + A₁⋅α⋅θ⋅λ - 2⋅A₁⋅α⋅θ + A₁⋅α⋅λ - 2⋅A₁⋅α - 2⋅R_F0
──────────────────────────────────────────────────────────────────────────────
                  R_M⋅θ⋅(λ - 2)
                                               A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
        2            2                       - ─────────── + ───────── + 2⋅A⋅α
⋅θ - R_M ⋅θ⋅λ + 2⋅R_M ⋅θ + R_M⋅R_W⋅θ⋅λ              2            2
──────────────────────────────────────, R_L: ─────────────────────────────────

                                                             A₁⋅α⋅θ⋅λ
⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅A⋅α - A₁⋅R_W⋅α⋅θ⋅λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α⋅θ
                                                                2
──────────────────────────────────────────────────────────────────────────────
                                     θ⋅(λ - 2)
   A₁⋅α⋅λ                   R_M⋅R_W⋅θ⋅λ⎫⎤
 + ────── - A₁⋅α - R_F0⋅θ + ───────────⎪⎥
     2                           2     ⎬⎥
───────────────────────────────────────⎪⎥
                                       ⎭⎦



"""




"""

;;;;;;;;;;;;;;;;;;;  ss5_7


ss5_7 = solve((eq1, eq2, eq3, eq4, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, P_B, R_L
                  , dict=True)



⎡⎧
⎢⎪
⎢⎨                      R_F0⋅θ                        -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ
⎢⎪P_B: ────────────────────────────────────────, R_D: ────────────────────────
⎣⎩     -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ


 - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ        A⋅α + 2⋅A₁⋅α
─────────────────────────────────────────────────────────, R_D1: ─────────────
            θ⋅(λ - 2)


- θ⋅(A⋅α + 2⋅A₁⋅α - R_M)       -A⋅R_W⋅α⋅θ⋅λ + A⋅R_W⋅α⋅λ - 2⋅A₁⋅R_W⋅α⋅θ⋅λ + 2⋅A
────────────────────────, R_E: ───────────────────────────────────────────────
     θ                                                               λ - 2


₁⋅R_W⋅α⋅λ - 2⋅R_F0⋅θ + R_M⋅R_W⋅θ⋅λ
──────────────────────────────────, R_E1: -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R

             A⋅R_W⋅α⋅θ⋅λ   A⋅R_W⋅α⋅λ
           - ─────────── + ───────── + 2⋅A⋅α⋅θ⋅λ - 4⋅A⋅α⋅θ - A⋅α⋅λ + 2⋅A⋅α - A
                  2            2
_M⋅θ, R_L: ───────────────────────────────────────────────────────────────────

                           A₁⋅α⋅θ⋅λ            A₁⋅α⋅λ                   R_M⋅R_
₁⋅R_W⋅α⋅θ⋅λ + A₁⋅R_W⋅α⋅λ + ──────── - A₁⋅α⋅θ + ────── - A₁⋅α - R_F0⋅θ + ──────
                              2                  2                           2
──────────────────────────────────────────────────────────────────────────────
   θ⋅(λ - 2)
W⋅θ⋅λ⎫⎤
─────⎪⎥
     ⎬⎥
─────⎪⎥
     ⎭⎦




"""



"""

;;;;;;;;;;;;;;;;;;;;;;;;;;; ss5_8

ss5_8 = solve(( eq2, eq3, eq4, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, R_L
                  , dict=True)


⎡⎧
⎢⎪
⎢⎨     -R_F0⋅(2⋅P_B - R_W⋅λ)                A⋅α            2⋅A₁⋅α
⎢⎪R_D: ──────────────────────, R_D1: -A⋅α + ─── - 2⋅A₁⋅α + ────── + R_M, R_E:
⎣⎩          P_B⋅(λ - 2)                      θ               θ


-R_F0⋅θ⋅(2⋅P_B - R_W⋅λ)
────────────────────────, R_E1: -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ, R_L:
      P_B⋅(λ - 2)
                                                       A₁⋅P_B⋅α⋅θ⋅λ
 2⋅A⋅P_B⋅α⋅θ⋅λ - 4⋅A⋅P_B⋅α⋅θ - A⋅P_B⋅α⋅λ + 2⋅A⋅P_B⋅α + ──────────── - A₁⋅P_B⋅α
                                                            2
 ─────────────────────────────────────────────────────────────────────────────
                                                            P_B⋅θ⋅(λ - 2)
     A₁⋅P_B⋅α⋅λ                           R_F0⋅R_W⋅θ⋅λ⎫⎤
⋅θ + ────────── - A₁⋅P_B⋅α - P_B⋅R_F0⋅θ + ────────────⎪⎥
         2                                     2      ⎬⎥
──────────────────────────────────────────────────────⎪⎥
                                                      ⎭⎦




"""



"""

;;;;;;;;;;;;;;;;;; ss5_9


ss5_9 = solve((eq2, eq3, eq4, beq2, beq3, beq4)
                  , R_D1, R_E1, R_D, R_E, R_L, R_F0
                  , dict=True)

[]




"""




"""

;;;;;;;;;;;;;;;;;;;;;;;; ss5_11





⎡⎧     -R_F0⋅(2⋅P_B - R_W⋅λ)         A⋅α + 2⋅A₁⋅α - θ⋅(A⋅α + 2⋅A₁⋅α - R_M)
⎢⎨R_D: ──────────────────────, R_D1: ─────────────────────────────────────, R_
⎣⎩          P_B⋅(λ - 2)                                θ
   -R_F0⋅θ⋅(2⋅P_B - R_W⋅λ)                                                 ⎫⎤
E: ────────────────────────, R_E1: -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ⎬⎥
         P_B⋅(λ - 2)                                                       ⎭⎦




"""



"""
ss5_12 = solve((eq2, eq3, beq2, beq3)
                   , R_D1, R_E1, R_E
                   , dict=True)



ss5_13 = solve((eq2, eq3, eq4, beq3)
                   , R_D1, R_E1, R_D, R_E
                   , dict=True)



⎡⎧     -R_F0⋅(2⋅P_B - R_W⋅λ)         A⋅α + 2⋅A₁⋅α - θ⋅(A⋅α + 2⋅A₁⋅α - R_M)
⎢⎨R_D: ──────────────────────, R_D1: ─────────────────────────────────────, R_
⎣⎩          P_B⋅(λ - 2)                                θ
   -R_F0⋅θ⋅(2⋅P_B - R_W⋅λ)                                                 ⎫⎤
E: ────────────────────────, R_E1: -A⋅α⋅θ + A⋅α - 2⋅A₁⋅α⋅θ + 2⋅A₁⋅α + R_M⋅θ⎬⎥
         P_B⋅(λ - 2)                                                       ⎭⎦



"""

