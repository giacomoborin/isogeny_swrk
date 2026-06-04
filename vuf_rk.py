import hashlib
import logging
from random import randint
import time

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)
from sage.all import *
proof.all(False)

import vuf.params as params
from vuf.vuf import VUF
import vuf.ec as ec
from rerandom_keys_vuf import H_path_vuf, H_path_sk_vuf, path_to_new_sk_vuf, gen_two_tors_basis
from update_sign_vuf import update_sig
import vuf.hd as hd

# TODO: non-adaptable version
class VUFRK(VUF):
    def __init__(self, adaptable=True):
        self.adaptable = adaptable

        # Basis of 2^f on E0
        self.B_E0_2f = (params.U0, params.V0)

        # Exponent params
        self.f1 = params.f1
        self.f = params.f
        self.ucmp = self.f - self.f1

        # Number of iterations
        # Length of walk
        l = log(params.p,2) + 255 +3*log(params.N, 2)

        r = ceil(l/self.ucmp)
        while r*self.ucmp < l + log(r*self.ucmp,2.):
            r += 1
        print(f'Set {r = }, {self.ucmp = }, for total walk length {ceil(l)}')

        self.r = r
        logger.info(f'Number of pushforwards: {self.r}, for total walk length {ceil(l)} and u={self.ucmp}')

        self.sk, self.pk = self.keygen()

        # Basis of 2^f on E_pk
        _, _, _, U1, V1 =self.sk
        self.B_Epk_2f = (U1, V1)
        self.B_Epk_2f_plain = gen_two_tors_basis(self.pk[0], params.f)
        self.M_Epk_2f = ec.ChangeOfBasis(self.B_Epk_2f, self.B_Epk_2f_plain)


    

    def get_randomness(self):
        rr = []
        for i in range(self.r):
            rr.append(randint(0,2**(self.ucmp)-1))
        return rr

    def expand(self,rr):
        E, R, S = self.pk
        return H_path_vuf(E, self.B_Epk_2f_plain, [R,S], rr, e=self.ucmp)

    def rand_pk(self,rr):
        """
        Online key rerandomization.
        """
        logger.info('Starting online public key update')
        _t0 = time.time()
        _, Es, _, pts = self.expand(rr)
        E_pk = Es[-1]
        R, S = pts[-1]
        pk = (E_pk, R, S)
        _t1 = time.time()
        logger.info(f'Total public key update time: {_t1-_t0:.3f}s')
        return pk
    
    def rand_keys(self,rr):
        """
        Offline key rerandomization.
        """
        
        logger.info('Starting secret key update')
        _t0 = time.time()

        E_pk, I_sk, M_sk, _, _ = self.sk
        _, R, S = self.pk

        pushed_pts, cKs, lms, Est, Bst, ee = H_path_sk_vuf(E_pk,self.B_Epk_2f_plain,(R,S),rr,self.ucmp)
    
        _t1 = time.time()
        logger.info(f'- Isogeny path (expand) done: {_t1-_t0:.3f}s')

        sk = path_to_new_sk_vuf(I_sk, M_sk,self.B_E0_2f, self.M_Epk_2f, cKs, lms, Est, Bst, ee)
        
        pk = (sk[0], pushed_pts[0], pushed_pts[1])


        _t2 = time.time()
        logger.info(f'- Secret key computation done: {_t2-_t1:.3f}s')

        return sk, pk

def IsogVerify_rk(P1, Q1, P3, Q3, v):

    K = ((P1, P3), (Q1, Q3))
    phi = hd.Dim2Iso(K, params.f1)

    if v[0].is_isomorphic(phi.codomain()[0]) and v[1].is_isomorphic(phi.codomain()[1]):
        return phi
    if v[1].is_isomorphic(phi.codomain()[0]) and v[0].is_isomorphic(phi.codomain()[1]):
        return phi
    return False


def VUFRK_verify(msg, sigma, pk, UV):
    """
    Verification
    Input:
    - msg: message
    - sigma = (pi, v): signature
    - pk: public key
    Output:
    - boolean verification
    """
    pi, v = sigma 
    E, R, S = pk
    if not R or not S or params.N * R or params.N * S:
        return False
    U3, V3 = pi

    # TODO: U and V should be a deterministic basis, we pull it from the secret
    # key for now
    U, V = UV

    f_inx = int(hashlib.shake_256(msg).hexdigest((params.f1 // 8) + 1), 16)
    r, s = 1, f_inx % params.N
    Rx = r*R + s*S
    if not Rx:
        return False

    # Verify the dim2 iso
    cof = 2**(params.f - params.f1)
    P1 = cof * U
    Q1 = cof * V

    phi = IsogVerify_rk(P1, Q1, U3, V3, v)
    if not phi:
        return False

    assert r != 0, 'r = 0, restart' # should not happen, we set r=1
    Ex, phiR, phiS = hd.EmbeddedIsogeny(phi, params.f1, params.N, (R, S))

    Qphi = r * phiR + s * phiS
    return Qphi == 0

def VUFRK_update(sigma, pk, B_Epk_2f, rr):
    r"""
    Signature update function
    Input:
    - sigma: signature
    - pk: public key
    - B_Epk_2f: full torsion basis of E_pk[2^f]
    - rr: randomness
    Output:
    - updated signature
    """
    logger.info('Starting signature update')
    _t0 = time.time()

    E_pk, R, S = pk
    # pi, v = sigma

    _t1 = time.time()

    u = params.f-params.f1

    Ks, Es, Bs, pushed_pts, Kds = H_path_vuf(E_pk, B_Epk_2f, [R,S], rr, e=u, dual=True)

    pk_rr = (Es[-1], pushed_pts[-1][0], pushed_pts[-1][1])

    _t2 = time.time()
    logger.info(f'- Isogeny path (expand) done: {_t2-_t1:.3f}s')

    (U3, V3), _ = sigma 

    # Reformating
    E3 = U3.curve()
    sigma = (
            (U3, V3), E3
        )

    # TODO: B_Epk_2f_rr should really be deterministically generated when used for verify
    #  but for now we output it
    B_Epk_2f_rr, sigma_rr = update_sig(E_pk,B_Epk_2f,sigma,params.N,rr, Es, Bs, Kds, u)

    _t3 = time.time()

    logger.info(f'- Pushforward done: {_t3-_t2:.3f}s')

    logger.info(f'Total signature update: {_t3-_t0:.3f}s')


    return pk_rr, sigma_rr, B_Epk_2f_rr


if __name__ == "__main__":
    # Add logging
    logger.setLevel(logging.DEBUG)

    alice = VUFRK()
    rr = alice.get_randomness()
    set_random_seed(rr[0])
    alice.rand_pk(rr)
    sk_new, pk_new = alice.rand_keys(rr)

    msg = b'Hello world'
    tstart = time.time()
    sigma = alice.VUFeval(msg)
    _t = time.time()
    logger.info(f'Signing time: {_t-tstart:.3f}s')

    tstart = time.time()
    out = alice.verify(sigma[0], msg, sigma[1])
    _t = time.time()
    logger.info(f'Verifying the original signature: {out}')
    logger.info(f'Verification time: {_t-tstart:.3f}s')


    pk_rr, sigma_rr, B_Epk_2f_rr = VUFRK_update(sigma, alice.pk, alice.B_Epk_2f, rr)


    state = VUFRK_verify(msg, sigma_rr, pk_rr, B_Epk_2f_rr)
    logger.info(f'\nVerification of updated signature: {state}')
