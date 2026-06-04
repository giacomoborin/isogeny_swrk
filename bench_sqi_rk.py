import time

from sage.all import randint

from sqisign_prism_v2.sqisign import SQIsign_verify
import sqisign_prism_v2.params as params
from sqisign_rk import SQIsignRK

if __name__ == "__main__":
    #assert False, 'run benchmarks with `sage --python -O bench_sqi.py [level]` to skip assertions'

    lvl = 1#int(sys.argv[1]) if len(sys.argv) > 1 else 1
    params.set_sqi_params(lvl)
    print(f'Running SQIsignRK lvl {lvl}')

    n_runs = 10

    tot_times = [0, 0, 0, 0, 0] # [keygen, sign, verification]

    for i in range(n_runs):
        print(f'\t- run {i+1}/{n_runs}')

        _t0 = time.time()
        sqi = SQIsignRK()
        _t1 = time.time()
        print(f'\t\t- Keygen time: {_t1-_t0}')

        rr = sqi.get_randomness()
        pk = sqi.rand_pk(rr)
        _t2 = time.time()
        print(f'\t\t- Online public key randomization time: {_t2-_t1}')

        sqi.rand_keys(rr)
        assert pk == sqi.pk
        _t3 = time.time()
        print(f'\t\t- Offline key randomization time: {_t3-_t2}')


        msg = f'Hello world {randint(1, 100)}'
        sigma = sqi.sign(msg)

        _t4 = time.time()
        print(f'\t\t- Signing time: {_t4-_t3}')

        out = SQIsign_verify(msg, sigma, sqi.pk)
        assert out
        _t5 = time.time()
        print(f'\t\t- Verification time: {_t5-_t4}')

        tot_times[0] += _t1 - _t0
        tot_times[1] += _t2 - _t1
        tot_times[2] += _t3 - _t2
        tot_times[3] += _t4 - _t3
        tot_times[4] += _t5 - _t4


    names = ['Keygen', 'Pk update', 'Sk update', 'Sign', 'Verify']
    print(f'Done')
    for i in range(5):
        tt = tot_times[i] / n_runs
        print(f'{names[i]}:\t{tt:.3f}s')
    tt = (sum(tot_times)) / n_runs
    print(f'Total time:\t{tt:.3f}s')

    print(f'Ratio Pk update / Keygen: {tot_times[1]/tot_times[0]:.2f}')
    print(f'Ratio Sk update / Keygen: {tot_times[2]/tot_times[0]:.2f}')

