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
import sqisign_prism_v2.ec as ec
from sqisign_prism_v2.sqisign import SQIsign, SQIsign_verify
from rerandom_keys import H_path, path_to_new_sk

class SQIsignRK(SQIsign):
    def get_randomness(self):
        # using ucmp = f, rrand = 2
        return (randint(0,2**(params.f)-1),randint(0,2**(params.f)-1))

    def expand(self,rr):
        return H_path(self.pk,rr,params.f)

    def rand_pk(self,rr):
        """
        Online key rerandomization.
        """
        logger.info('Starting online public key update')
        _t0 = time.time()
        _, Es, Bs = self.expand(rr)
        pk = (Es[1], Bs[1])
        _t1 = time.time()
        logger.info(f'Total public key update time: {_t1-_t0:.3f}s')
        return pk

    def rand_keys(self,rr=None):
        """
        Offline key rerandomization.
        """
        logger.info('Starting offline public and secret key update')
        _t0 = time.time()
        if rr is None:
            logger.info('Needs to genrate randomness')
            rr = self.get_randomness()
        _, Es, Bs = self.expand(rr)
        _t1 = time.time()
        logger.info(f'- Isogeny path (expand) done: {_t1-_t0:.3f}s')
        (P_pk_bis, Q_pk_bis) = ec.TorsionBasis(Es[-1])
        self.pk = (Es[-1], (P_pk_bis, Q_pk_bis))
        _t2 = time.time()
        logger.info(f'- Public key (torsion basis) computation done: {_t2-_t1:.3f}s')
        self.sk = path_to_new_sk(self.sk,rr,Es,Bs,params.f)
        _t3 = time.time()
        logger.info(f'- Secret key computation done: {_t3-_t2:.3f}s')
        logger.info(f'Total public and secret key update time: {_t3-_t0:.3f}s')

if __name__=="__main__":
    # Add logging
    logger.setLevel(logging.DEBUG)

    lvl = 1
    params.set_sqi_params(lvl)
    alice = SQIsignRK()

    alice.rand_keys()

    msg = 'Hello world'
    sigma = alice.sign(msg)

    out = SQIsign_verify(msg, sigma, alice.pk)
    print(f'Verification: {out}')






