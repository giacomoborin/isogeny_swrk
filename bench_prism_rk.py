import time
import sys

from sage.all import randint

import sqisign_prism_v2.params as params
from prism_rk import PRISMRK, PRISMRK_update, PRISMRK_verify

if __name__ == "__main__":
    #assert False, 'run benchmarks with `sage --python -O bench_sqi.py [level]` to skip assertions'

    lvl = 5#int(sys.argv[1]) if len(sys.argv) > 1 else 1
    params.set_prism_params(lvl)
    print(f'Running PRISMRK lvl {lvl}')

    n_runs = 10

    tot_times = [0, 0, 0, 0, 0, 0]

    for i in range(n_runs):
        print(f'\t- run {i+1}/{n_runs}')

        _t0 = time.time()
        prism = PRISMRK()
        _t1 = time.time()
        print(f'\t\t- Keygen time: {_t1-_t0}')

        rr = prism.get_randomness()
        pk, basis_plain = prism.rand_pk(rr)
        _t2 = time.time()
        print(f'\t\t- Online public key randomization time: {_t2-_t1}')

        msg = f'Hello world {randint(1, 100)}'
        sigma = prism.sign(msg)
        _t3 = time.time()
        print(f'\t\t- Signing time: {_t3-_t2}')

        sigma_rr = PRISMRK_update(msg,sigma,prism.pk,prism.B_pk_plain,rr)
        _t4 = time.time()
        print(f'\t\t- Signature update time: {_t4-_t3}')

        out = PRISMRK_verify(msg, sigma_rr, pk)
        assert out
        _t5 = time.time()
        print(f'\t\t- Verification time: {_t5-_t4}')

        prism.rand_keys(rr)
        assert pk == prism.pk
        assert basis_plain == prism.B_pk_plain
        _t6 = time.time()
        print(f'\t\t- Offline key randomization time: {_t6-_t5}')

        tot_times[0] += _t1 - _t0
        tot_times[1] += _t2 - _t1
        tot_times[2] += _t3 - _t2
        tot_times[3] += _t4 - _t3
        tot_times[4] += _t5 - _t4
        tot_times[5] += _t6 - _t5


    names = ['Keygen', 'Pk update', 'Sign', 'Adapt', 'Verify', 'Sk update']
    print(f'Done')
    for i in range(6):
        tt = tot_times[i] / n_runs
        print(f'{names[i]}:\t{tt:.3f}s')
    tt = (sum(tot_times)) / n_runs
    print(f'Total time:\t{tt:.3f}s')

    print(f'Ratio Pk update / Keygen: {tot_times[1]/tot_times[0]:.2f}')
    print(f'Ratio Sk update / Keygen: {tot_times[5]/tot_times[0]:.2f}')
    print(f'Ratio Adapt / Sign: {tot_times[3]/tot_times[2]:.2f}')