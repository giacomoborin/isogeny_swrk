from sage.all import *
import hashlib
import logging
import time 
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

try:
    from . import params, hd
    from . import quaternions as qt
    from . import qlapoti as qlpt
    from .theta.theta_structures.couple_point import CouplePoint
except ImportError:
    from theta.theta_structures.couple_point import CouplePoint
    import params, hd
    import quaternions as qt
    import qlapoti as qlpt

def IsogenyRepresentation(I, J, E, U, V, U0, V0):
    """
    Algorithm 5
    """
    nK = 2**params.f1 - params.N
    K = qt.RandomIdealGivenNorm(nK, False)
    L = K.intersection(I.intersection(J))

    _, U3, V3 = qlpt.IdealToIsogeny(L, U0, V0)

    s = inverse_mod(I.norm(), 2**params.f)
    U3, V3 = s*U3, s*V3

    cof = 2**(params.f - params.f1)
    K = (
            (cof*params.N*U, cof*U3),
            (cof*params.N*V, cof*V3)
        )
    phi = hd.Dim2Iso(K, params.f1)
    # TODO: detect which is Ex and which Fx
    return phi.codomain(), U3, V3

def IsogVerify(P1, Q1, P3, Q3, v):
    K = ((params.N*P1, P3), (params.N*Q1, Q3))
    phi = hd.Dim2Iso(K, params.f1)

    if v[0].is_isomorphic(phi.codomain()[0]) and v[1].is_isomorphic(phi.codomain()[1]):
        return phi
    if v[1].is_isomorphic(phi.codomain()[0]) and v[0].is_isomorphic(phi.codomain()[1]):
        return phi
    return False

class VUF():
    def __init__(self):
        sk, pk = self.keygen()
        self.sk = sk
        self.pk = pk

    def keygen(self):
        _t0 = time.time()
        Nsk = next_prime(randint(1, params.p**2))
        I = qt.RandomIdealGivenNorm(Nsk, True)

        # Evaluate the ideal on the twist
        _, U1, V1, P, Q = qlpt.IdealToIsogeny(
                I, params.U0, params.V0, params.P0, params.Q0)
        logger.info('Twisted evaluation done')

        # Scaling N-torsion
        a, b, c, d = [randint(1, params.N) for _ in range(4)]
        R, S = a * P + b * Q, c * P + d * Q
        M = [a, b, c, d]
        s = inverse_mod(Nsk, 2**params.f)
        U1, V1 = s*U1, s*V1

        Epk = P.curve()
        _t1 = time.time()
        logger.info(f'Key gen took: {_t1-_t0:.3f}s')
        # TODO: scale U0, V0 from a deterministic basis
        return (Epk, I, M, U1, V1), (Epk, R, S)

    def VUFeval(self, x):
        Epk, I, M, U1, V1 = self.sk
        a, b, c, d = M

        # TODO: better hashing
        f_inx = int(hashlib.shake_256(x).hexdigest((params.f1 // 8) + 1), 16)
        r, s = 1, f_inx % params.N

        gamma = (r*a + s*c) + (r*b + s*d)*params.kappa
        Ix = params.O0 * params.alpha.conjugate() * gamma.conjugate() + params.O0 * params.N
        Ex, U3, V3 = IsogenyRepresentation(
                I, Ix, Epk, U1, V1, params.U0, params.V0)

        # Sanity
        e1 = U1.weil_pairing(V1, 2**params.f)
        e3 = U3.weil_pairing(V3, 2**params.f)
        ee = e1 ** (params.N * (2**params.f1 - params.N))
        assert e3 == ee

        # TODO: Ex is now (Ex, Fx)
        pi = (U3, V3)
        v = Ex
        logger.info('Evaluation done')
        return pi, v

    def verify(self, pi, x, v):

        E, R, S = self.pk   
        if not R or not S or params.N * R or params.N * S:
            return False
        U3, V3 = pi

        # TODO: U and V should be a deterministic basis, we pull it from the secret
        # key for now
        U, V = self.sk[3:]

        f_inx = int(hashlib.shake_256(x).hexdigest((params.f1 // 8) + 1), 16)
        r, s = 1, f_inx % params.N
        Rx = r*R + s*S
        if not Rx:
            return False

        # Verify the dim2 iso
        cof = 2**(params.f - params.f1)
        P1 = cof * U
        Q1 = cof * V
        P3 = cof * U3
        Q3 = cof * V3
        phi = IsogVerify(P1, Q1, P3, Q3, v)
        if not phi:
            return False

        assert r != 0, 'r = 0, restart' # should not happen, we set r=1
        Ex, phiR, phiS = hd.EmbeddedIsogeny(phi, params.f1, params.N, (R, S))

        Qphi = r * phiR + s * phiS
        return Qphi == 0

if __name__ == "__main__":
    # Add logging
    logger.setLevel(logging.DEBUG)

    rr = randint(1, 2**64)

    set_random_seed(rr)

    logger.info(f'Done with precomputation. Running with seed {rr}')
    vuf = VUF()

    msg = b'Hello world'
    pi, v = vuf.VUFeval(msg)

    state = vuf.verify(pi, msg, v)
    print(f'\nVerification: {state}')

