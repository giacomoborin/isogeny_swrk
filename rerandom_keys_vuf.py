from sage.all import *
proof.all(False)

import vuf.params as params
import vuf.quaternions as qt
import vuf.ec as ec
import vuf.qlapoti as qlpt

# from montgomery_isogenies.isogenies_x_only import (isogeny_from_kernel_x_only, 
# 	evaluate_isogeny_x_only, 
# 	evaluate_isogeny_x_only_one_point)


def exp_order_two(P):
	r"""
	Returns e such that order(P)=2**e.
	"""

	e = 0
	Q = P
	while Q!=0:
		Q = 2*Q
		e += 1

	return e

def gen_two_tors_basis(E, e):
	r"""
	Very naive function to generate a basis of E[2^e]. 
	"""
	cof = (params.p+1)*(params.p-1) // 2**e
	while True:
		P = E.random_point() * cof
		if P * 2**(e-1):
			break

	while True:
		Q = E.random_point() * cof
		if Q * 2**(e-1):
			e2 = P.weil_pairing(Q, 2**e)
			if e2**(2**(e-1)) != 1:
				break

	return P, Q

def H_path_vuf(E_pk,B_Epk_2f,pts,rr,e=params.f, dual=False):
    r"""
    Outputs r isogenies of degree 2^e, using randomness rr of length r.
    Adapted to be used in VUFRK

    Input:
    - E_pk: the public key curve 
    - B_Epk_2f: a basis of E_pk[2^f].
    - pts: a list of points on the public key curve.
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
    - pushed_pts: images of the input points under the isogeny path.
    - Kds (if dual): list of phi_i(Q_i) generating \ker(\hat{phi_i}) 
    (0 <= i <= r-1).
    """

    r = len(rr)

    Ks = []

    Edom = E_pk
    Es = []
    P, Q = B_Epk_2f
    P = 2**(params.f-e)*P
    Q = 2**(params.f-e)*Q

    Bs = []
    Kds = []
    pushed_pts = [pts]
    for i in range(r):
        K = P + rr[i]*Q
        # Replaced because it was too slow
        K._order = ZZ(2**e)
        phi = Edom.isogeny(K, model="montgomery")
        Ecod = phi.codomain()
        # phi, Ecod = isogeny_from_kernel_x_only(Edom, ZZ(2**e), K)
        Es += [Ecod]
        pushed_pts += [[phi(pt) for pt in pushed_pts[-1]]]
        if dual:
            # Replaced because it was too slow
            Kds +=[phi(Q)]
            # Kds +=[evaluate_isogeny_x_only_one_point(phi,Q)]
        P, Q = gen_two_tors_basis(Ecod, params.f)
        Bs += [(P,Q)]
        P = 2**(params.f-e)*P
        Q = 2**(params.f-e)*Q

        Edom = Ecod

    # todo add mult for pushed_pts
    if dual:
        return Ks, Es, Bs, pushed_pts, Kds
    else:
        return Ks, Es, Bs, pushed_pts


def KernelDecomposedToIdealLT_vuf(c1, c2, e=params.f):
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
    M = a*params.mat_1+(-1-b)*params.mat_i+2*b*params.mat_ij2+b*params.mat_1k2
    X = ec.EvalMatrix(M)
    P = c1*X[0] + c2*X[1]
    assert 2**(params.f-e)*P == 0

    ker = params.B((a + b/2, -1, b, b/2))
    I = params.O0 * ker + params.O0 * (2**e)
    assert I.norm() == 2**e
    return I


def H_path_sk_vuf(E_pk, B_Epk_2f, pts_pk, rr, e):

	r"""Analogue of H_path adapted to the secret key update in VURFRK.
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
	- B_Epk_2f: basis of E_pk[2^f].
	- pts_pk: list of points on E_pk to be pushed through the isogeny path.
	- rr: a list of r random integers between 0 and 2^e-1, determining the
	public key rerandomizing 2-isogeny path.
	- e: integer smaller or equal to f = v_2(p+1).

	Output:
	- pushed_pts: list of the input points pushed through the isogeny path.
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
	P, Q = B_Epk_2f
	P = 2**(params.f-e)*P
	Q = 2**(params.f-e)*Q
    
	Es = []
	Bs = []
	pushed_pts = pts_pk

	for i in range(r):
		K = P + rr[i]*Q
		# Replaced because it was too slow
		K._order = ZZ(2**e)
		phi = Edom.isogeny(K, model="montgomery")
		Ecod = phi.codomain()
		# phi, Ecod = isogeny_from_kernel_x_only(Edom, ZZ(2**e), K)
		Es += [Ecod]
		L_phis.append(phi)
		P, Q = gen_two_tors_basis(Ecod, params.f)
		Bs += [(P,Q)]
		P *= 2**(params.f-e)
		Q *= 2**(params.f-e)
		Edom = Ecod
		# Push pk points through
		pushed_pts = [phi(pt) for pt in pushed_pts]

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
			P, Q = B_Epk_2f
		else:
			P, Q = Bs[i*n_main-1]
		P *= 2**(params.f-e*n_main)
		Q *= 2**(params.f-e*n_main)
		
		phiP = P
		phiQ = Q
		for j in range(n_main):
			# Replaced because it was too slow
			phiP = L_phis[i*n_main+j](phiP)
			phiQ = L_phis[i*n_main+j](phiQ)
			# phiP, phiQ = evaluate_isogeny_x_only(L_phis[i*n_main+j], phiP, phiQ)
		
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
		phiP = L_phis[r_main*n_main+j](phiP)
		phiQ = L_phis[r_main*n_main+j](phiQ)
		# phiP, phiQ = evaluate_isogeny_x_only(L_phis[r_main*n_main+j], phiP, phiQ)


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

	return pushed_pts, cKs, lms, Est, Bst, ee


def path_to_new_sk_vuf(I_sk, M_sk, B2, M2, cKs, lms, Est, Bst, ee):
	r"""
	Updates VUF secret key from a precomputed random 2-isogeny
	path (with H_path_sk_vuf, see above).

	Input:
	- I_sk: I_sk being the secret key ideal corresponding to phi_sk: E_0 --> E_pk
	- M_sk: change of basis matrix for R,S in E_pk[N], part of secret key
	- B2: a basis (U1, V1) of E_pk[2^f].
	- M2: change of basis matrix from (phi_sk(U_0), phi_sk(V_0)) to 
	(U1, V1), where (U_0, V_0) is a basis of E_0[2^f] and 
	(U1, V1) is a basis of E_pk[2^f].
	- cKs, lms, Est, Bst, ee: outputs of H_path_sk_prism, see above.

	Output:
	- E: new public key.
	- J: new secret key ideal.
	- M_sk: updated secret change of basis matrix for images of R,S in E[N]
	- (P, Q): deterministic basis of E[2^f].
	"""
	r = len(ee)
	I = I_sk
	M = Matrix(Zmod(2**ee[0]), M2)
	M_sk_new = Matrix(Zmod(params.N), 2, M_sk)

	for i in range(r):
		E = Est[i]
		(P, Q) = Bst[i]
		# Computing I
		NI = ZZ(I.norm())
		a, b, c, d = M.list()
		# If M = [[a,b], [c,d]] and K = [cKs[i][0]]*P + [cKs[i][1]]*Q
		# then \hat{\phi_{J}}(K) = [N]([cKs[i][0]*a+cKs[i][1]*c]P0+[cKs[i][0]*b+cKs[i][1]*d]Q0)
		c1 = NI*(cKs[i][0]*a+cKs[i][1]*c)
		c2 = NI*(cKs[i][0]*b+cKs[i][1]*d)


		pullback_IK = KernelDecomposedToIdealLT_vuf(c1, c2, ee[i]-2*lms[i])
		IK = ZZ(2**lms[i])*qt.Pushforward(pullback_IK, I)
		IpIK = I*IK

		## Computing J and phi_J
		J, gamma  = qt.RandomEquivalentPrimeIdeal(IpIK, iso=True)
		U, V = B2
		EJ, UJ, VJ = qlpt.IdealToIsogeny(J, U, V)
		assert EJ.j_invariant() == E.j_invariant()

		eta = EJ.isomorphism_to(E)
		PJ, QJ = eta(UJ), eta(VJ)

		M = ec.ChangeOfBasis((PJ, QJ), (P, Q))
		I = J

		# Update M_sk according to gamma
		c1,c2,c3,c4 = qt.DecomposeQuaternionAlongBasis(gamma, [params.B(1), params._i, params.zeta, params._i*params.zeta], params.N)
		M_gamma = c1*params.mat_1_N + c2*params.mat_iota_N + c3*params.mat_zeta_N + c4*params.mat_iotazeta_N
		NI_inv = inverse_mod(NI, params.N)
		M_sk_new = NI_inv * M_sk_new * M_gamma
		
	# Reformat M_sk_new as a list of 4 integers
	M_sk_new = [ZZ(x) for x in M_sk_new.list()]

	return E, J, M_sk_new, P, Q 

