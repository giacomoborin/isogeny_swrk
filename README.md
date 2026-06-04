# Isogeny signatures with Randomizable keys in Sage

SageMath implementation of isogeny-based signatures with rerandomizable keys
(SQIsign-RK, PRISM-RK, DeuringVUF-RK), for the paper:

> "Isogeny-based Signatures with Randomizable Keys"
> Andrea Basso, Giacomo Borin, Maria Corte-Real Santos, Pierrick Dartois, Riccardo Invernizzi, Luciano Maino, Robi Pedersen and Michel Seck

**Disclaimer**: The provided code is a proof of concept implementation, neither optimized for performance nor for security.

**Authors**:
- Riccardo Invernizzi
- Maria Corte-Real Santos

### Dependencies

- SQIsign and PRISM original implementation were obtained from [KULeuven-COSIC/sqisign_prism_v2](https://github.com/KULeuven-COSIC/sqisign_prism_v2), included in the submodule `sqisign_prism_v2`.
- Since SageMath native functions were deemed too slow, one dimensional 2-isogeny chain computations use Giacomo Pope's x-only arithmetic codebase (in `montgomery_isogenies` and `utilities`), available at [GiacomoPope/KummerIsogeny](https://github.com/GiacomoPope/KummerIsogeny).


## Usage - SQIsign-RK

To test SQIsign original (`KeyGen`, `Sign`, `Verif`) and new (`RandPK`, `RandPK`) procedures at once, open a `sage` terminal and type:
```python
load("sqisign_rk.py")
```
This should execute a single run at NIST level `lvl=1` (NIST levels 3 or 5 may be set
by changing the `lvl` variable in `sqisign_rk.py`).

To benchmark 10 full executions of SQIsign-RK (e.g. to try and reproduce timings from Table 3 of the paper), type:
```python
load("bench_sqi_rk.py")
```
By default, the benchmark is for NIST level `lvl=1` but the `lvl` variable may be changed accordingly in `bench_sqi_rk.py` to benchmark other levels.

## Usage - PRISM-RK

To test PRISM original (`KeyGen`, `Sign`, `Verif`) and new (`RandPK`, `RandPK`, `Adapt`) procedures at once, open a `sage` terminal and type:
```python
load("prism_rk.py")
```
This should execute a single run at NIST level `lvl=1` (NIST levels 3 or 5 may be set
by changing the `lvl` variable in `sqisign_rk.py`).

To benchmark 10 full executions of PRISM-RK (e.g. to try and reproduce timings from Table 3 of the paper), type:
```python
load("bench_prism_rk.py")
```
By default, the benchmark is for NIST level `lvl=1` but the `lvl` variable may be changed accordingly in `bench_prism_rk.py` to benchmark other levels.

## Usage - VUF-RK

In the folder `vuf` we provide a SageMath/Python implementation of the DeuringVUF signature. It can be run on its own by opening a `sage` terminal in that folder and typing:
```python
load("vuf.py")
```

To test DeruingVUF original (`KeyGen`, `Sign`, `Verif`) and new (`RandPK`, `RandPK`, `Adapt`) procedures at once, open a `sage` terminal and type:
```python
load("vuf_rk.py")
```
This should execute a single run at NIST level `lvl=1`.

To benchmark 10 full executions of VUF-RK (e.g. to try and reproduce timings from Table 3 of the paper), type:
```python
load("bench_vuf_rk.py")
```

## Project structure

### Main code

- `sqisign_rk.py`: main implementation of SQIsign-RK inheriting from `sqisign_prism_v2/sqisign.py`
- `prism_rk.py`: main implementation of PRISM `sqisign_prism_v2/prism.py`
- `vuf_rk.py`: main implementation of VUF-RK inheriting from `vuf/vuf.py`
- `rerandom_keys.py`: functions to rerandomize keys (`RandPK`, `RandSK`)
- `update_sign.py`: functions to update signatures (`Adapt`)
- `rerandom_keys_vuf.py`: functions to rerandomize keys (`RandPK`, `RandSK`) adapted to VUF-RK (as arithmetic is over Fp4).
- `update_sign_vuf.py`: functions to update signatures (`Adapt`)  adapted to VUF-RK (as arithmetic is over Fp4).

### Libraries and functions

- `sqisign_prism_v2`: submodule imported from [KULeuven-COSIC/sqisign_prism_v2](https://github.com/KULeuven-COSIC/sqisign_prism_v2) and containing original SQIsign and PRISM implementations
- `montgomery_isogenies`, `utilities`: Giacomo Pope's x-only arithmetic codebase extracted from [GiacomoPope/KummerIsogeny](https://github.com/GiacomoPope/KummerIsogeny)

### Benchmarking

- `bench_sqi_rk.py`: benchmarking for SQIsign-RK
- `bench_prism_rk.py`: benchmarking for PRISM-RK
- `bench_vuf_rk.py`: benchmarking for VUF-RK
