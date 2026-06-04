from sage.all import *
proof.all(False)

#import sys
#sys.path.append("sqisign_prism_v2/")

import sqisign_prism_v2.params as params
import sqisign_prism_v2.quaternions as qt
import sqisign_prism_v2.ec as ec
import sqisign_prism_v2.qlapoti as qlpt
from sqisign_prism_v2.sqisign import SQIsign
from montgomery_isogenies.isogenies_x_only import (isogeny_from_kernel_x_only, 
	evaluate_isogeny_x_only, 
	evaluate_isogeny_x_only_one_point)



def H_path(pk, rr, ucmp=params.f, dual=False):
    r"""
    Outputs r isogenies of degree 2^e, using randomness rr of length r.

    Input:
    - pk = (E_pk, (P_pk, Q_pk)): public key in the SQIsign format, E_pk 
    being the public key curve and (P_pk, Q_pk) being a deterministic basis
    of E_pk[2^e].
    - rr: a list of r random integers between 0 and 2^e-1, determining the
    public key rerandomizing 2-isogeny path.
    - e: integer smaller or equal to f = v_2(p+1) (default value f).
    - dual: boolean value for optional output (Kds).

    Output:
    - Ks: list of kernels K_i:=P_i+[rr[i]]*Q_i, where (P_i, Q_i) is a
    deterministic basis of E_i[2^e] (0 <= i <= r-1, with E_0=E_pk).
    - Es: list of codomains of isogenies phi_i: E_i -> E_{i+1} with 
    kernel K_i = Ks[i] (0 <= i <= r-1).
    - Bs: list of deterministic basis (P_{i+1}, Q_{i+1}) of E_{i+1}[2^e]
    (0 <= i <= r-1).
    - Kds (if dual): list of phi_i(Q_i) generating \ker(\hat{phi_i}) 
    (0 <= i <= r-1).
    """
    E_pk, (P_pk, Q_pk) = pk

    r = len(rr)

    Ks = []

    Edom = E_pk
    P,Q = P_pk, Q_pk
    Es = []
    Bs = []
    Kds = []

    for i in range(r):
        K = P + rr[i]*Q
        # Replaced because it was too slow
        # K._order = ZZ(2**e)
        # phi = Edom.isogeny(K, model="montgomery")
        # Ecod = phi.codomain()
        phi, Ecod = isogeny_from_kernel_x_only(Edom, ZZ(2**ucmp), K)
        Es += [Ecod]
        if dual:
            # Replaced because it was too slow
            # Kds +=[phi(Q)]
            Kds +=[evaluate_isogeny_x_only_one_point(phi,Q)]
        P, Q = ec.TorsionBasis(Ecod)
        Bs += [(P,Q)]
        P = 2**(params.f-ucmp)*P
        Q = 2**(params.f-ucmp)*Q
        Edom = Ecod

    if dual:
        return Ks, Es, Bs, Kds
    else:
        return Ks, Es, Bs

def KernelDecomposedToIdealLT(c1, c2, e=params.f):
    """
    Given c1, c2 defining a point K = c1*P0 + c2*Q0 in E0[2^e] return the ideal
    corresponding to the isogeny with kernel K [Algorithm 3.17]

    Input:
    - c1, c2: integers defining a point K = c1*P0 + c2*Q0 in E0[2^e]
    - e: integer <= f (default f)

    Output:
    - I: a left O0-ideal corresponding to the isogeny of kernel K
    """
    v = vector(Zmod(2**params.f), (c1, c2))

    # theta = j + (1+k)/2 (and j = -i + 2*((i+j)/2))
    M_theta = -params.mat_i + 2*params.mat_ij2 + params.mat_1k2
    d1, d2 = M_theta.transpose() * v # c1*theta(P0) + c2*theta(Q0) = theta(P)

    M = Matrix(Zmod(2**params.f), 2, [c1, d1, c2, d2])

    a, b = M**-1 * params.mat_i.transpose() * v
    a, b = ZZ(a), ZZ(b)

    # Now aP + b*thetaP == eta*P
    # Sanity:
    import sqisign_prism_v2.ec as ec
    M = a*params.mat_1+(-1-b)*params.mat_i+2*b*params.mat_ij2+b*params.mat_1k2
    X = ec.EvalMatrix(M)
    P = c1*X[0] + c2*X[1]
    assert 2**(params.f-e)*P == 0

    ker = params.B((a + b/2, -1, b, b/2))
    I = params.O0 * ker + params.O0 * (2**e)
    assert I.norm() == 2**e
    return I


def path_to_new_sk(sk,rr,Es,Bs,e=params.f):
    r"""
    Updates SQIsign secret key from a precomputed random 2-isogeny
    path (with H_path, see above).

    Input:
    - sk = (E_pk, (P_pk, Q_pk), I_sk, M_sk): secret key in SQIsign format,
    E_pk being the public key curve, (P_pk, Q_pk) being a deterministic basis
    of E_pk[2^e], I_sk being the secret key ideal and M_sk being the change of
    basis matrix from (phi_sk(P_0), phi_sk(Q_0)) to (P_pk, Q_pk), where (P_0, Q_0)
    is the public deterministic basis of E_0[2^f] (f can be != e).
    - rr: a list of r random integers between 0 and 2^e-1, determining the
    public key rerandomizing 2-isogeny path.
    - Es (output of H_path): list of codomains of isogenies phi_i: E_i -> E_{i+1} 
    with kernel K_i = P_i+[rr[i]]*Q_i (0 <= i <= r-1).
    - Bs (output of H_path): list of deterministic basis (P_{i+1}, Q_{i+1}) 
    of E_{i+1}[2^e] (0 <= i <= r-1).
    - e: integer smaller or equal to f = v_2(p+1) (default value f).

    Output:
    - sk_new: new secret key after rerandomization in SQIsin format 
    (see sk above).
    """

    r = len(rr)
    _, _, I_sk, M_sk = sk
    I = I_sk
    M = Matrix(Zmod(2**e), M_sk)
    for i in range(r):
        E = Es[i]
        (P, Q) = Bs[i]
        # Computing I
        N = ZZ(I.norm())
        a, b, c, d = M.list()
        # If M = [[a,b], [c,d]] and K = P + [rr[i]]*Q
        # then \hat{\phi_{J}}(K) = [N]([a+rr[i]*c]P+[b+rr[i]*d]Q)
        c1 = N*(a+rr[i]*c)
        c2 = N*(b+rr[i]*d)

        pullback_IK = KernelDecomposedToIdealLT(c1, c2, e)
        # I.IK = I.([I]_*[I]^*IK) = I inter [I]^*IK 
        # = nrd(I)[I]^*IK + nrd([I]^*IK)I.
        I_IK = N*pullback_IK + ZZ(2**e)*I
        IK = qt.Pushforward(pullback_IK, I)
        IpIK = I*IK
        assert I_IK == IpIK


        ## Computing J and phi_J
        J = qt.RandomEquivalentPrimeIdeal(I_IK)

        assert J.is_left_equivalent(I_IK)

        EJ, PJ, QJ = qlpt.IdealToIsogeny(J)

        assert EJ.j_invariant() == E.j_invariant()

        eta = EJ.isomorphism_to(E)
        PJ, QJ = eta(PJ), eta(QJ)

        M = ec.ChangeOfBasis((PJ, QJ), (P, Q))
        I = J
        # assert ec.EvalMatrix(M, (PJ, QJ)) == (P, Q)
        # a, b, c, d = M.list() # M = [[a,b], [c,d]]


    sk_new = E, (P, Q), J, M

    return sk_new

def exp_order_two(P):
    r"""Returns e such that order(P)=2**e.
    """

    e = 0
    Q = P
    while Q!=0:
        Q = 2*Q
        e += 1

    return e

def H_path_sk_prism(E_pk,B_pk_plain,rr,e):
    r"""Analogue of H_path adapted to the secret key update of PRISM-RK.
    Computes r isogenies of degree 2^e, using randomness rr of length r
    and regroups them in r_main chunks psi_i of length n_main:

    E_{i*n_main-1} --phi_{i*n_main}--> E_{i*n_main} ... 
    E_{i*n_main+j-1} --phi_{i*n_main+j}--> E_{i*n_main+j} ...
    E_{(i+1)*n_main-2} --phi_{(i+1)*n_main-1}--> E_{(i+1)*n_main-1}

    and a last chunck psi_{r_main} of length n_last:

    E_{r_main*n_main-1} --phi_{r_main*n_main}--> E_{r_main*n_main} ... 
    E_{r_main*n_main+j-1} --phi_{r_main*n_main+j}--> E_{r_main*n_main+j} ...
    E_{r_main*n_main+n_last-2} --phi_{r_main*n_main+n_last-1}--> E_{r_main*n_main+n_last-1}

    so that n_main is the biggest integer such that n_main*e <= f in order to maximize the
    use of 2^f-torsion and minimize the number of ideal to isogeny translations in 
    path_to_new_sk_prism.

    Input:
    - E_pk: public key curve.
    - B_pk_plain = (P_pk_plain, Q_pk_plain): full torsion deterministic basis
    of E_pk[2^f] (f = v_2(p+1)).
    - rr: a list of r random integers between 0 and 2^e-1, determining the
    public key rerandomizing 2-isogeny path.
    - e: integer smaller or equal to f = v_2(p+1).

    Output:
    - cKs: list of r_main+1 integers a_i, b_i such that 
    ker(psi'_i)=[2^m_i]<[a_i]*P_i+[b_i]*Q_i>,
    where (P_i,Q_i) is a deterministic basis of the deg(psi_i)-torsion
    basis of the domain of psi_i, psi'_i is cyclic and psi_i = [2^m_i]psi'_i.
    - lms: list of the integers m_i from above.
    - Est: list of codomains of the psi_i.
    - Bst: list of deterministic basis of the 2^f-torsion of the codomain of psi_i.
    - ee: list of e_i such that deg(phi_i)=2^e_i
    """
    r = len(rr)

    L_phis = []

    Edom = E_pk
    P,Q = B_pk_plain
    Es = []
    Bs = []

    P *= 2**(params.f-e)
    Q *= 2**(params.f-e)

    for i in range(r):
        K = P + rr[i]*Q
        # Replaced because it was too slow
        # K._order = ZZ(2**e)
        # phi = Edom.isogeny(K, model="montgomery")
        # Ecod = phi.codomain()
        phi, Ecod = isogeny_from_kernel_x_only(Edom, ZZ(2**e), K)
        Es += [Ecod]
        L_phis.append(phi)
        P, Q = ec.TorsionBasis(Ecod)
        Bs += [(P,Q)]
        P = 2**(params.f-e)*P
        Q = 2**(params.f-e)*Q
        Edom = Ecod

    n_main = params.f//e
    r_main = r//n_main
    n_last = r%n_main

    ee = [n_main*e]*r_main
    ee.append(n_last*e)

    cKs = []
    lms = []
    Est = []
    Bst = []

    # Main iteration
    for i in range(r_main):
        if i == 0:
            P, Q = B_pk_plain
        else:
            P, Q = Bs[i*n_main-1]
        P *= 2**(params.f-e*n_main)
        Q *= 2**(params.f-e*n_main)

        phiP = P
        phiQ = Q
        for j in range(n_main):
            # Replaced because it was too slow
            # phiP = L_phis[i*n_main+j](phiP)
            # phiQ = L_phis[i*n_main+j](phiQ)
            phiP, phiQ = evaluate_isogeny_x_only(L_phis[i*n_main+j], phiP, phiQ)

        # ker(phi)=<a*P+b*Q>, a*phiP + b*phiQ = 0
        eP = exp_order_two(phiP)
        eQ = exp_order_two(phiQ)
        em = max(eP,eQ)
        dm = e*n_main-em
        cphiP = 2**(dm)*phiP
        cphiQ = 2**(dm)*phiQ
        if eP == em:
            a = ZZ(pari.elllog(Es[i*n_main+n_main-1],-cphiQ,cphiP,ZZ(2**(em))))
            b = ZZ(1)
        else:
            a = ZZ(1)
            b = ZZ(pari.elllog(Es[i*n_main+n_main-1],-cphiP,cphiQ,ZZ(2**(em))))

        assert a*cphiP + b*cphiQ == 0

        cKs.append((a,b))
        lms.append(e*n_main-em)
        Est.append(Es[i*n_main+n_main-1])
        Bst.append(Bs[i*n_main+n_main-1])

    # Last chunk of the chain
    P, Q = Bs[r_main*n_main-1]
    P *= 2**(params.f-e*n_last)
    Q *= 2**(params.f-e*n_last)

    phiP = P
    phiQ = Q
    for j in range(n_last):
        # Replaced because it was too slow
        # phiP = L_phis[r_main*n_main+j](phiP)
        # phiQ = L_phis[r_main*n_main+j](phiQ)
        phiP, phiQ = evaluate_isogeny_x_only(L_phis[r_main*n_main+j], phiP, phiQ)


    # ker(phi)=<a*P+b*Q>, a*phiP + b*phiQ = 0
    eP = exp_order_two(phiP)
    eQ = exp_order_two(phiQ)
    em = max(eP,eQ)
    dm = e*n_last-em
    cphiP = 2**(dm)*phiP
    cphiQ = 2**(dm)*phiQ
    if eP == em:
        a = ZZ(pari.elllog(Es[-1],-cphiQ,cphiP,ZZ(2**(em))))
        b = ZZ(1)
    else:
        a = ZZ(1)
        b = ZZ(pari.elllog(Es[-1],-cphiP,cphiQ,ZZ(2**(em))))

    assert a*cphiP + b*cphiQ == 0

    cKs.append((a,b))
    lms.append(e*n_last-em)
    Est.append(Es[-1])
    Bst.append(Bs[-1])

    return cKs, lms, Est, Bst, ee

def path_to_new_sk_prism(I_sk, M_sk_plain, cKs, lms, Est, Bst, ee):
    r"""
    Updates PRISM secret key from a precomputed random 2-isogeny
    path (with H_path_sk_prism, see above).

    Input:
    - I_sk: I_sk being the secret key ideal corresponding to phi_sk: E_0 --> E_pk
    - M_sk_plain: change of basis matrix from (phi_sk(P_0), phi_sk(Q_0)) to 
    (P_pk, Q_pk), where (P_0, Q_0) is the public deterministic basis of 
    E_0[2^f] and (P_pk, Q_pk) is a deterministic basis of E_pk[2^f].
    - cKs, lms, Est, Bst, ee: outputs of H_path_sk_prism, see above.

    Output:
    - E: new public key.
    - (P, Q): deterministic basis of E[2^f].
    - J: new secret key ideal.
    - M: change of basis matrix from (phi_J(P_0), phi_J(Q_0)) to (P, Q).
    """
    r = len(ee)
    I = I_sk
    M = Matrix(Zmod(2**ee[0]), M_sk_plain)
    for i in range(r):
        E = Est[i]
        (P, Q) = Bst[i]
        # Computing I
        N = ZZ(I.norm())
        a, b, c, d = M.list()
        # If M = [[a,b], [c,d]] and K = [cKs[i][0]]*P + [cKs[i][1]]*Q
        # then \hat{\phi_{J}}(K) = [N]([cKs[i][0]*a+cKs[i][1]*c]P0+[cKs[i][0]*b+cKs[i][1]*d]Q0)
        c1 = N*(cKs[i][0]*a+cKs[i][1]*c)
        c2 = N*(cKs[i][0]*b+cKs[i][1]*d)


        pullback_IK = KernelDecomposedToIdealLT(c1, c2, ee[i]-2*lms[i])
        IK = ZZ(2**lms[i])*qt.Pushforward(pullback_IK, I)
        IpIK = I*IK

        ## Computing J and phi_J
        J = qt.RandomEquivalentPrimeIdeal(IpIK)

        EJ, PJ, QJ = qlpt.IdealToIsogeny(J)

        assert EJ.j_invariant() == E.j_invariant()

        eta = EJ.isomorphism_to(E)
        PJ, QJ = eta(PJ), eta(QJ)

        M = ec.ChangeOfBasis((PJ, QJ), (P, Q))
        I = J
        # assert ec.EvalMatrix(M, (PJ, QJ)) == (P, Q)
        # a, b, c, d = M.list() # M = [[a,b], [c,d]]


    sk_new = E, (P, Q), J, M

    return sk_new


if __name__=="__main__":
    lvl = 1
    params.set_sqi_params(lvl)
    alice = SQIsign()

    rr = (randint(0,2**(params.f)-1),randint(0,2**(params.f)-1))

    Ks, Es, Bs = H_path(alice.pk,rr)
    sk2 = path_to_new_sk(alice.sk,rr,Es,Bs)

