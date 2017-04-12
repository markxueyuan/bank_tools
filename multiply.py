from sympy.solvers.solveset import nonlinsolve
from sympy import Symbol, symbols, init_printing,\
    simplify, integrate, latex, solve_poly_system,\
    Function, Derivative, apart, cancel, expand, collect,\
    Integral

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


# symbols as undefined function

f = symbols("f", cls=Function)
r = symbols('r', cls=Function)

##################### Equations #######################

# Firm, FOC
feq1 = P - R_L/r(L)


# Mutual Fund, FOC

meq1 = R_B - R_F0


# HH, FOC

eq1 = R_D1 * θ * P_B - (R_F0 * θ)
eq2 = R_E1 - R_D1 * θ
eq3 = -2*R_F0*θ + 2*R_E - λ*R_D*θ + R_F0 * θ * λ * R_W / P_B
eq4 = -2*R_F0*θ + R_D*θ + (1-λ)*R_D*θ + R_F0 * θ * λ * R_W / P_B

## Bank, FOC

beq1 = R_I - R_M - P*f(I)
beq2 = -2*R_L + R_M**2 + R_M*R_I
beq3 = R_M - R_D1(1-α*(A+A1)+A1*α)-R_E1*(α*(A+A1)+A1*α)/R_E1


e0 = 2*α*A
e1 = (1+2*λ*R_W)*α*A

beq4 = 2*R_L - 2*R_E*e0/R_E - R_D*(1-e0/R_E)- R_D1*(A1*α/R_E1) \
       - R_E1*(A1*α/R_E1)




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
        ⎛  λ    ⎞   R_D1⋅R_W⋅λ
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





if __name__ == "__main__":
    main()














