import hashlib
import logging
import time

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)
from sage.all import *
proof.all(False)

import sqisign_prism_v2.params as params
import sqisign_prism_v2.quaternions as qt
import sqisign_prism_v2.ec as ec
import sqisign_prism_v2.qlapoti as qlpt
import sqisign_prism_v2.hd as hd
from sqisign_prism_v2.prism import PRISM, PRISM_verify
from rerandom_keys import H_path, H_path_sk_prism, path_to_new_sk_prism
from update_sign import update_sig

def hash_to_prime_msg(msg, r=None):
    """
    Hash (message || r) into an (a-1)-bit number
    randomly generating r until the output is prime.
    """

    if type(msg) == str: msg = msg.encode()
    if b'&&' in msg:
        # TODO: better domain separation
        raise ValueError('invalid message')

    msg = msg + b'&&'

    if r:
        chl = msg + r
        h = hashlib.sha256(chl).digest()
        q = int.from_bytes(h, 'big') % (2**params.a)
        if is_pseudoprime(q):
            return q, r
        raise ValueError("invalid salt provided")

    while True:
        r = randint(0, 2**params.rb)
        r = int(r).to_bytes((params.rb+7)// 8, 'big')
        # cnt = str(counter).encode()
        chl = msg + r
        h = hashlib.sha256(chl).digest()
        q = int.from_bytes(h, 'big') % (2**params.a)
        if is_pseudoprime(q):
            return q, r

class PRISMRK(PRISM):
    def __init__(self):
        # Additional params
        self.ucmp = params.f-params.a
        self.r = ceil(2*log(params.p)/(log(2)*self.ucmp))
        # Key generation
        self.pk, self.sk, self.M_sk_plain, self.B_pk_plain = self.key_gen()

    def key_gen(self):
        """
        Key generation
        Output:
        - sk: secret key
        - pk: public key
        """
        logger.info('Starting keygen')
        _t0 = time.time()

        I_sk = qt.RandomIdealGivenNorm(params.D_mix, True)
        assert I_sk.left_order() == params.O0 and I_sk.norm() == params.D_mix

        I_sk = qt.RandomEquivalentPrimeIdeal(I_sk)

        _t1 = time.time()
        logger.info(f'- Ideal generation done: {_t1-_t0:.3f}s')

        E_pk, phi_P0, phi_Q0 = qlpt.IdealToIsogeny(I_sk)

        _t2 = time.time()
        logger.info(f'- Qlapoti done: {_t2-_t1:.3f}s')

        P_pk, Q_pk = ec.TorsionBasis(E_pk)
        # We do not want to multiply by 2**(params.f - params.a) for now, because 
        # the change of basis matrix from (phi_P0, phi_Q0) to (P_pk, Q_pk) (with full
        # torsion) will be needed for update
        #P_pk *= 2**(params.f - params.a)
        #Q_pk *= 2**(params.f - params.a)
        M_sk_plain = ec.ChangeOfBasis((phi_P0, phi_Q0), (P_pk, Q_pk))
        assert ec.EvalMatrix(M_sk_plain, (phi_P0, phi_Q0)) == (P_pk, Q_pk)

        basis_plain = (P_pk, Q_pk)

        P_pk *= 2**(params.f - params.a)
        Q_pk *= 2**(params.f - params.a)
        M_sk = matrix(Integers(2**params.f),
            [[M_sk_plain[0,0]*2**(params.f - params.a),
            M_sk_plain[0,1]*2**(params.f - params.a)],
            [M_sk_plain[1,0]*2**(params.f - params.a),
            M_sk_plain[1,1]*2**(params.f - params.a)]])

        pair_pk = pari.ellweilpairing(E_pk, P_pk, Q_pk, 2**params.a)
        pk = (E_pk, (P_pk, Q_pk), pair_pk)
        sk = (E_pk, I_sk, M_sk)

        assert basis_plain[0]*2**params.a != 0

        _t3 = time.time()
        logger.info(f'- Change of basis done: {_t3-_t2:.3f}s')
        logger.info(f'Keygen done: {_t3-_t0:.3f}s')

        return pk, sk, M_sk_plain, basis_plain

    def get_randomness(self):
        rr = []
        for _ in range(self.r):
            rr.append(randint(0,2**(self.ucmp)-1))
        return rr

    def expand(self,rr,dual=False):
        E_pk, B_pk, _ = self.pk
        P = 2**(params.a-self.ucmp)*B_pk[0]
        Q = 2**(params.a-self.ucmp)*B_pk[1]

        return H_path((E_pk, (P,Q)),rr,self.ucmp,dual)

    def rand_pk(self,rr):
        """
        Online key rerandomization.
        """
        logger.info('Starting public key update')

        _t0 = time.time()
        _, Es, Bs = self.expand(rr)
        _t1 = time.time()
        logger.info(f'- Isogeny path (expand) done: {_t1-_t0:.3f}s')

        E_pk_bis = Es[-1]
        B_pk_plain_bis = Bs[-1]

        P_pk_bis = 2**(params.f-params.a)*B_pk_plain_bis[0]
        Q_pk_bis = 2**(params.f-params.a)*B_pk_plain_bis[1]

        pair_pk_bis = pari.ellweilpairing(E_pk_bis, P_pk_bis, Q_pk_bis, 2**params.a)
        _t2 = time.time()
        logger.info(f'- Weil pairing done: {_t2-_t1:.3f}s')

        pk = (E_pk_bis, (P_pk_bis,Q_pk_bis), pair_pk_bis)
        logger.info(f'Total public key update time: {_t2-_t0:.3f}s')

        return pk, B_pk_plain_bis

    def rand_keys(self,rr):
        """
        Offline key rerandomization.
        """
        logger.info('Starting secret key update')

        _t0 = time.time()
        E_pk, I_sk, _ = self.sk

        cKs, lms, Est, Bst, ee = H_path_sk_prism(E_pk,self.B_pk_plain,rr,self.ucmp)
        _t1 = time.time()
        logger.info(f'- Isogeny path (expand) done: {_t1-_t0:.3f}s')

        E_pk_bis, (P_pk_plain_bis, Q_pk_plain_bis), I_sk_bis, M_sk_plain_bis = path_to_new_sk_prism(I_sk, self.M_sk_plain, cKs, lms, Est, Bst, ee)

        self.M_sk_plain = M_sk_plain_bis
        self.B_pk_plain = (P_pk_plain_bis, Q_pk_plain_bis)

        M_sk_bis = matrix(Integers(2**params.f),
            [[M_sk_plain_bis[0,0]*2**(self.ucmp),
            M_sk_plain_bis[0,1]*2**(self.ucmp)],
            [M_sk_plain_bis[1,0]*2**(self.ucmp),
            M_sk_plain_bis[1,1]*2**(self.ucmp)]])

        self.sk = (E_pk_bis, I_sk_bis, M_sk_bis)

        _t2 = time.time()
        logger.info(f'- Secret key computation done: {_t2-_t1:.3f}s')

        P_pk_bis = 2**(self.ucmp)*P_pk_plain_bis
        Q_pk_bis = 2**(self.ucmp)*Q_pk_plain_bis

        pair_pk_bis = pari.ellweilpairing(E_pk_bis, P_pk_bis, Q_pk_bis, 2**params.a)

        self.pk = (E_pk_bis, (P_pk_bis,Q_pk_bis), pair_pk_bis)

        _t3 = time.time()
        logger.info(f'- Public key computation done: {_t3-_t2:.3f}s')
        logger.info(f'Total secret key update time: {_t3-_t0:.3f}s')

    def sign(self, msg):
        """
        Signing.
        Input:
        - msg: the message
        Output:
        - sig: a valid signature
        """
        logger.info('Starting signing')
        _t0 = time.time()

        E_pk, I_sk, M_sk = self.sk

        # Hash the message to a prime
        q, r = hash_to_prime_msg(msg)
        assert q < 2**params.a

        _t1 = time.time()
        logger.info(f'- Hashing done: {_t1-_t0:.3f}s')

        # Construct the response to the challenge
        n_rsp = q*(2**params.a - q)
        I_rsp = qt.RandomIdealGivenNorm(n_rsp, False)
        I_rsp = qt.Pushforward(I_rsp, I_sk)
        I_cra = I_sk * I_rsp

        _t2 = time.time()
        logger.info(f'- Response ideal done: {_t2-_t1:.3f}s')

        # Compute the corresponding isogeny
        E_rsp, P_cra, Q_cra = qlpt.IdealToIsogeny(I_cra)

        # Scale torsion and multiplication by q^-1
        q_inv = inverse_mod(q, 2**params.a)
        P_rsp, Q_rsp = ec.EvalMatrix(q_inv*M_sk, basis=(P_cra, Q_cra))

        pts_rsp = (P_rsp, Q_rsp)
        sigma = (
            E_rsp, pts_rsp, r
        )

        _t3 = time.time()
        logger.info(f'- Response isogeny done: {_t3-_t2:.3f}s')
        logger.info(f'Signature done: {_t3-_t0:.3f}s')
        return sigma

def PRISMRK_verify(msg, sigma, pk):
    """
    Verification
    Input:
    - msg: message
    - sigma: signature
    - pk: public key
    Output:
    - boolean verification
    """
    logger.info('Starting verification')
    _t0 = time.time()

    # (P_pk, Q_pk) are 2^a-torsion
    E_pk, (P_pk, Q_pk), pair_pk = pk
    E_rsp, pts_rsp, r = sigma
    q, r1 = hash_to_prime_msg(msg, r = r)
    assert r1 == r # Hashing gives a prime

    P_rsp, Q_rsp = pts_rsp

    # Compute a 2D isogeny to check the response
    K = ((P_pk, P_rsp), (Q_pk, Q_rsp))
    Phi = hd.Dim2Iso(K, params.a)

    _t1 = time.time()
    logger.info(f'- 2D isogeny done: {_t1-_t0:.3f}s')

    # Check the degree using pairings (à la SQIsign2D-East)
    P = hd.CouplePoint(P_pk, E_rsp(0))
    Q = hd.CouplePoint(Q_pk, E_rsp(0))

    P1, Q1 = Phi(P)[0], Phi(Q)[0]

    pair = pari.ellweilpairing(P1.curve(), P1, Q1, 2**params.a)
    pair_q = pair_pk ** q
    pair_qinv = pair_q ** (-1)

    if pair in [pair_q, pair_qinv]:
        _t2 = time.time()
        logger.info(f'- Degree checking done: {_t2-_t1:.3f}s')
        logger.info(f'Verification done: {_t2-_t0:.3f}s')
        return True
    logger.error('Degree not matching')
    return False

def PRISMRK_update(msg,sigma,pk,B_pk_plain,rr):
    r"""
    Signature update function
    Input:
    - msg: message
    - sigma: signature
    - pk: public key
    - B_pk_plain: full torsion basis of E_pk[2^f]
    - rr: randomness
    Output:
    - updated signature
    """
    logger.info('Starting signature update')
    _t0 = time.time()

    E_pk, (P_pk, Q_pk), _ = pk
    _, _, r = sigma

    q, r1 = hash_to_prime_msg(msg, r = r)
    _t1 = time.time()

    u = params.f-params.a

    Ks, Es, Bs, Kds = H_path((E_pk,(2**(params.a-u)*P_pk,2**(params.a-u)*Q_pk)), rr,ucmp=u, dual=True)
    _t2 = time.time()
    logger.info(f'- Isogeny path (expand) done: {_t2-_t1:.3f}s')

    sigma_rr = update_sig(pk,B_pk_plain,sigma,q,rr, Es, Bs, Kds, u)
    _t3 = time.time()
    logger.info(f'- Pushforward done: {_t3-_t2:.3f}s')

    logger.info(f'Total signature update: {_t3-_t0:.3f}s')

    return sigma_rr






if __name__ == "__main__":
    # Add logging
    logger.setLevel(logging.DEBUG)

    # Setting parameters
    lvl = 1#int(sys.argv[1]) if len(sys.argv) > 1 else 1
    params.set_prism_params(lvl)

    # We need to add u to parameters
    alice = PRISMRK()

    rr = alice.get_randomness()
    pk, basis_plain = alice.rand_pk(rr)

    msg = 'Hello world'
    sigma = alice.sign(msg)
    sigma_rr = PRISMRK_update(msg,sigma,alice.pk,alice.B_pk_plain,rr)

    out = PRISMRK_verify(msg, sigma_rr, pk)
    print(f'Verification: {out}')

    alice.rand_keys(rr)

    assert pk == alice.pk
    assert basis_plain == alice.B_pk_plain

    msg = 'Hello world'
    sigma = alice.sign(msg)

    out = PRISMRK_verify(msg, sigma, alice.pk)
    print(f'Verification: {out}')
    
