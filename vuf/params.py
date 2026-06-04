from sage.all import (
        next_prime, ceil, log, sqrt, ZZ, QuaternionAlgebra, GF, EllipticCurve,
        Matrix, Zmod, floor, next_prime, randint, Integers, gcd, kronecker
)

try:
    from . import ec
except ImportError:
    import ec

try:
    from . import quaternions
except ImportError:
    import quaternions

# New set of parameters
# print('Starting setup')
p = 611 * 2**262 - 1
N = 255893511804611991997964586367829068761994890719409435889324626467241377
f = 262
f1 = 238 # 2^f1 > N

# f = 500
# p = 27 * 2**f - 1
# N = 502812446935834963588720280717315223454589901895482372644296761715125242633595193
# ell = 2
# f1 = ZZ(N).bit_length()



# Fields and curves
Fp = GF(p)
Fp2, Fp2_i = GF(p**2, name="i", modulus=[1, 0, 1]).objgen()
Fp4, Fp4_t = Fp2.extension(2, name="t").objgen()

E0 = EllipticCurve(Fp4, [1, 0])
assert E0.order() == (p+1)**2*(p-1)**2
# E0 = EllipticCurve(Fp2, [1, 0])
# E0.set_order((p + 1) ** 2)

B = QuaternionAlgebra(-1, -p)
_i, _j, _k = B.gens()
O0 = B.maximal_order(order_basis=(B(1), _i, (_i+_j)/2, (1-_k)/2))

# Quadratic twist
# d_tw = 1 + 2*Fp2_i
# E0t = EllipticCurve(Fp2, [d_tw**2, 0])
# E0t.set_order((1 - p) ** 2)

iota = lambda P: P.curve()(-P[0], P[1] * Fp2_i)
frob = lambda P: P.curve()(P[0]**p, P[1]**p)
# k -> iota(frob(P))
# om = (1 + 2*Fp2_i)**((1-p) // 2)
# iota_t = lambda P: P.curve()(-P[0], P[1] * Fp2_i)
# frob_t = lambda P: P.curve()(om**2 * P[0]**p, om**3 * P[1]**p)

def eval_end(alpha, P, N=None):
    """
    Evaluate alpha (assumed to be integral) on P on the twist of E0. Optionally
    the order of P can be provided.
    """
    if not N:
        N = P.order()
    a, b, c, d = [x % N for x in alpha]
    # TODO: precompute the action of the endomorphisms
    return a * P + b * iota(P) + c * frob(P) + d * iota(frob(P))

# def eval_end_twist(alpha, P, N=None):
#     """
#     Evaluate alpha (assumed to be integral) on P on the twist of E0. Optionally
#     the order of P can be provided.
#     """
#     if not N:
#         N = P.order()
#     a, b, c, d = [x % N for x in alpha]
#     # TODO: precompute the action of the endomorphisms
#     return a * P + b * iota_t(P) + c * frob_t(P) + d * iota_t(frob_t(P))

# Sanity
# for _ in range(10):
#     P = E0t.random_point()
#     assert iota_t(P) in E0t and frob_t(P) in E0t
#     assert -p * P == frob_t(frob_t(P))

# Basis of the 2-torsion on E0
cof = (p+1)*(p-1) // 2**f
while True:
    U0 = E0.random_point() * cof
    if U0 * 2**(f-1):
        break

while True:
    V0 = E0.random_point() * cof
    if V0 * 2**(f-1):
        e2 = U0.weil_pairing(V0, 2**f)
        if e2**(2**(f-1)) != 1:
            break

# print('2-torsion basis found')

# Precomputing endomosphisms actions on 2 torsion
# i : (x, y) -> (-x, iy)
iota = lambda P: P.curve()(-P[0], P[1] * Fp2_i)

# pi : (x, y) -> (x^p, y^p)
frob = lambda P: P.curve()(P[0]**p, P[1]**p)

# Compute the points
U0_i, V0_i = iota(U0), iota(V0)

Fp8 = Fp4.extension(2, name="s")
E0_t = E0.change_ring(Fp8)
U0_t, V0_t = E0_t(U0), E0_t(V0)

U0_2 = U0_t.division_points(2)[0]
V0_2 = V0_t.division_points(2)[0]

U0_ij2 = E0(iota(U0_2) + frob(U0_2))
V0_ij2 = E0(iota(V0_2) + frob(V0_2))

U0_1k2 = E0(U0_2 + iota(frob(U0_2)))
V0_1k2 = E0(V0_2 + iota(frob(V0_2)))

# Compute matrices
mat_1 = Matrix(Zmod(2**f), 2, [1, 0, 0, 1])
mat_i = ec.ChangeOfBasis((U0, V0), (U0_i, V0_i))
mat_ij2 = ec.ChangeOfBasis((U0, V0), (U0_ij2, V0_ij2))
mat_1k2 = ec.ChangeOfBasis((U0, V0), (U0_1k2, V0_1k2))


# # For signature:

# # Generate alpha such that gcd(n(alpha), N^2) = N
# # print('Finding endomorphisms')
# while True:
#     b, c, d = [randint(1, N) for _ in range(3)]
#     n_res = b**2 + p * (c**2 + d**2)
#     if not Integers(N)(-n_res).is_square():
#         continue
#     alpha = B([Integers(N)(-n_res).sqrt(), b, c, d])
#     assert alpha.reduced_norm() % N == 0
#     if gcd(alpha.reduced_norm(), N**2) == N:
#         break


# # print(f'{alpha = }')
# IP0 = O0 * alpha.conjugate() + O0 * N

# # Compute kappa
# while True:
#     kappa = B([randint(1, N) for _ in range(4)])
#     if gcd(kappa.reduced_norm(), N) != 1:
#         continue
#     if kappa - list(kappa)[0] in IP0:
#         continue
#     break
# cof = (p-1)*(p+1)//N
# while True:
#     R0 = E0.random_point() * cof
#     if not R0:
#         continue # N is prime
#     P0 = eval_end(alpha, R0, N)
#     # P0 = eval_end_twist(alpha, R0, N)
#     if not P0:
#         continue
#     assert eval_end(alpha.conjugate(), P0, N) == 0
#     break

# For randomised signature
kappa = _i
n = 0
ZZN = Zmod(N)
while True:
    if kronecker(-p-n**2, N) == 1:
        lmbda = ZZN(-p-n**2).sqrt()
        lmbda = ZZ(lmbda)
        # print(f'{n = }, {lmbda = }')
        break
    n += 1
zeta = _j + n*_i
alpha = zeta+lmbda
# print(f'{alpha = }')

cof = (p-1)*(p+1)//N
while True:
    R0 = E0.random_point() * cof
    if not R0:
        continue # N is prime
    P0 = eval_end(alpha, R0, N)
    # P0 = eval_end_twist(alpha, R0, N)
    if not P0:
        continue
    assert eval_end(alpha.conjugate(), P0, N) == 0
    break

IP0 = O0 * (_j - lmbda) + O0 * N

# print(f'{P0 = }')
Q0 = eval_end(kappa, P0, N)
# Q0 = eval_end_twist(kappa, P0, N)
# print(f'{Q0 = }')

m1,m2,m3,m4 = quaternions.DecomposeQuaternionAlongBasis(zeta*_i, [B(1), _i, zeta, _i*zeta], N)
mat_1_N = Matrix(Zmod(N), 2, [1,0,0,1])
mat_iota_N = Matrix(Zmod(N), 2, [0,1,-1,0])
mat_zeta_N = Matrix(Zmod(N), 2, [lmbda, 0, m1+lmbda*m3, m2+lmbda*m4])
mat_iotazeta_N = mat_iota_N*mat_zeta_N

# Sanity

assert P0 * N == Q0 * N == 0
eN0 = P0.weil_pairing(Q0, N)
assert eN0 != 1 and eN0**N == 1

# Stuff from the sqisign routines
QUAT_prime_cofactor = next_prime(2**ceil(log(p, 2)))
QUAT_equiv_bound_coeff = 64

# print('Done')
