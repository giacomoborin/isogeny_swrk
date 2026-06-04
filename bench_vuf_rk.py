import time
import sys

from sage.all import randint

import vuf.params as params

from vuf_rk import VUFRK, VUFRK_update, VUFRK_verify

if __name__ == "__main__":
    #assert False, 'run benchmarks with `sage --python -O bench_sqi.py [level]` to skip assertions'

    n_runs = 10

    tot_times = [0, 0, 0, 0, 0, 0]

    for i in range(n_runs):
        print(f'\t- run {i+1}/{n_runs}')

        _t0 = time.time()
        vuf = VUFRK()
        _t1 = time.time()
        print(f'\t\t- Keygen time: {_t1-_t0}')

        rr = vuf.get_randomness()
        pk = vuf.rand_pk(rr)
        _t2 = time.time()
        print(f'\t\t- Online public key randomization time: {_t2-_t1}')

        msg = b'Hello world'
        sigma = vuf.VUFeval(msg)
        _t3 = time.time()
        print(f'\t\t- Signing time: {_t3-_t2}')

        pk_rr, sigma_rr, B_Epk_2f_rr = VUFRK_update(sigma, vuf.pk, vuf.B_Epk_2f, rr)
        _t4 = time.time()
        print(f'\t\t- Signature update time: {_t4-_t3}')

        pi, v = sigma_rr
        out = VUFRK_verify(pk_rr, pi, msg, v, B_Epk_2f_rr)
        assert out
        _t5 = time.time()
        print(f'\t\t- Verification time: {_t5-_t4}')

        vuf.rand_keys(rr)
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
