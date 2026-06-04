from sage.all import *
proof.all(False)

from time import time

import vuf.params as params
import vuf.quaternions as qt
import vuf.ec as ec
import vuf.qlapoti as qlpt
import vuf.hd as hd
from rerandom_keys_vuf import H_path_vuf, gen_two_tors_basis

def evaluate_basis(E_start,E_sigma,B_start_plain,B_start,B_sigma,q):
	r"""
	Input:
	- E_start: domain Montgomery curve
	- E_sigma: codomain Montgomery curve
	- B_start_plain: deterministic basis of E_start[2**f]
	- B_start = [2**(f-a)]*B_start_plain
	- B_sigma = [1/q]*sigma(B_start), with sigma: E_start -> E_sigma
	- q: integer such that deg(sigma) = q(2^a-q)
	Output:
	- sigma(B_start_plain)
	"""

	P_start_plain, Q_start_plain = B_start_plain
	P_start, Q_start = B_start
	P_sigma, Q_sigma = B_sigma

	# Compute a 2D isogeny of kernel {(P,[1/q]*sigma(P))| P\in E_start[2**a]}
	# sigma = phi2*phi1 and deg(phi1) = q with the notations 
	# from [https://eprint.iacr.org/2025/135.pdf, Theorem 1]
	Ker = ((P_start, P_sigma), (Q_start, Q_sigma))
	#assert P_start.order() == 2**params.a
	#assert Q_start.order() == 2**params.a
	#print(factor(P_sigma.order()))
	#assert P_sigma.order() == 2**params.a
	#assert Q_sigma.order() == 2**params.a
	#print(pari.ellweilpairing(E_start, P_start, Q_start, 2**params.a)*pari.ellweilpairing(E_sigma, P_sigma, Q_sigma, 2**params.a))
	#print(pari.ellweilpairing(E_start, P_start, Q_start, 2**params.a)/pari.ellweilpairing(E_sigma, P_sigma, Q_sigma, 2**params.a))
	#assert pari.ellweilpairing(E_start, P_start, Q_start, 2**params.a)*pari.ellweilpairing(E_sigma, P_sigma, Q_sigma, 2**params.a) == 1
	Phi = hd.Dim2Iso(Ker, params.f1)

	# Evaluating phi1(B_start_plain) = (phi1_P, phi1_Q)
	phi1_P, mpsi2_P = Phi(hd.CouplePoint(P_start_plain, E_sigma(0)))
	phi1_Q, mpsi2_Q = Phi(hd.CouplePoint(Q_start_plain, E_sigma(0)))

	# Evaluating \hat{phi2}(B_sigma_plain) = (hphi2_P, hphi2_Q)

	P_sigma_plain, Q_sigma_plain = gen_two_tors_basis(E_sigma, params.f)
	# ec.TorsionBasis(E_sigma)

	hphi2_P, psi1_P = Phi(hd.CouplePoint(E_start(0),P_sigma_plain))
	hphi2_Q, psi1_Q = Phi(hd.CouplePoint(E_start(0),Q_sigma_plain))

	# Correct basis order with pairings
	pair_start = pari.ellweilpairing(E_start, P_start_plain, Q_start_plain, 2**params.f)
	pair_phi1 =  pari.ellweilpairing(phi1_P.curve(), phi1_P, phi1_Q, 2**params.f)

	pair_start_q = pair_start**q
	pair_start_q_inv = 1/pair_start_q

	if pair_phi1 not in [pair_start_q, pair_start_q_inv]:
		#print("Swap")
		phi1_P = mpsi2_P
		phi1_Q = mpsi2_Q
		hphi2_P = psi1_P
		hphi2_Q = psi1_Q

		pair_phi1 =  pari.ellweilpairing(phi1_P.curve(), phi1_P, phi1_Q, 2**params.f)
	
		pair_start_q = pair_start**q
		pair_start_q_inv = 1/pair_start_q

	# Correct basis sign with pairings
	if pair_phi1 == pair_start_q_inv:
		#print("Minus phi1_Q")
		phi1_Q = -phi1_Q
	elif pair_phi1 != pair_start_q:
		raise ValueError("Wrong input degree q.")

	pair_sigma = pari.ellweilpairing(E_sigma, P_sigma_plain, Q_sigma_plain, 2**params.f)
	pair_hphi2 = pari.ellweilpairing(hphi2_P.curve(), hphi2_P, hphi2_Q, 2**params.f)

	pair_sigma_2amq = pair_sigma**(2**params.f1-q)
	pair_sigma_2amq_inv = 1/pair_sigma_2amq
	
	if pair_hphi2 == pair_sigma_2amq_inv:
		#print("Minus hphi2_Q")
		hphi2_Q = -hphi2_Q
	elif pair_hphi2 != pair_sigma_2amq:
		raise ValueError("Wrong input degree 2**a-q.")

	# Composition phi2*phi1(B_start)
	M_dlog = ec.ChangeOfBasis((hphi2_P,hphi2_Q),(phi1_P,phi1_Q))

	a,b,c,d = M_dlog.list()
	sigma_P = (2**params.f1-q)*(a*P_sigma_plain+b*Q_sigma_plain)
	sigma_Q = (2**params.f1-q)*(c*P_sigma_plain+d*Q_sigma_plain)

	pair_sigma = pari.ellweilpairing(E_sigma, sigma_P, sigma_Q, 2**params.f)
	
	#assert pair_sigma == pair_start**(q*(2**params.a-q))

	return sigma_P, sigma_Q

def push_forward(E_start,E_sigma,B_start_plain,B_start,B_sigma,q,E_push,B_push_plain,r,K_dual,e_phi):
	r"""
	Computes the pushforward of a 2D-representation (as in Algorithm 2 of the paper).

	Input:
	- E_start: domain Montgomery curve
	- E_sigma: codomain Montgomery curve
	- B_start_plain = (P_start_plain, Q_start_plain): deterministic basis of E_start[2**f]
	- B_start = [2**(f-a)]*B_start_plain
	- B_sigma = [1/q]*sigma(B_start), with sigma: E_start -> E_sigma
	- q: integer such that deg(sigma) = q(2^a-q)
	- E_push: codomain of an isogeny phi: E_start-> E_push
	- r: integer such that ker(phi)=<[2**(f-e_phi)]*(P_start_plain+r*Q_start_plain)>
	- K_dual: kernel generator of \hat{phi} ([2**(f-e_phi)]*phi(Q_start_plain))
	- e_phi: integer such that deg(phi)=2**e_phi
	Output:
	- [1/q][phi]_*(sigma)(B_push)
	"""
	
	P_push, Q_push = B_push_plain

	# Replaced because it was too slow
	K_dual._order = ZZ(2**e_phi)
	phi_dual = E_push.isogeny(K_dual, model="montgomery")
	eta = phi_dual.codomain().isomorphism_to(E_start)
	hphi_P, hphi_Q = eta(phi_dual(P_push)), eta(phi_dual(Q_push))
	# phi_dual, E_start_bis = isogeny_from_kernel_x_only(E_push, ZZ(2**e_phi), K_dual)
	# eta = E_start_bis.isomorphism_to(E_start)
	# hphi_P, hphi_Q = evaluate_isogeny_x_only(phi_dual, P_push, Q_push)
	# hphi_P, hphi_Q = eta(hphi_P), eta(hphi_Q)
	

	# DLOG matrix of phi_dual(B_push) in B_start_plain
	M_dlog = ec.ChangeOfBasis(B_start_plain,(hphi_P, hphi_Q))
	# phi_dual(P_push) = [a]*P_start_plain + [b]*P_start_plain
	# phi_dual(Q_push) = [c]*P_start_plain + [d]*P_start_plain
	a,b,c,d = M_dlog.list()

	## Compute sigma(B_start_plain)
	sigma_P,sigma_Q = evaluate_basis(E_start,E_sigma,B_start_plain,B_start,B_sigma,q)
	#print(f'sigma_P.order() = {factor(sigma_P.order())}')
	#print(f'sigma_Q.order() = {factor(sigma_Q.order())}')

	## Compute [sigma]_*(phi)
	# Image of ker(phi) via sigma
	sigma_K = 2**(params.f-e_phi)*(sigma_P+r*sigma_Q)
	
	# Replaced because it was too slow
	sigma_K._order = ZZ(2**e_phi)
	push_phi = E_sigma.isogeny(sigma_K, model="montgomery")
	# push_phi, _ = isogeny_from_kernel_x_only(E_sigma, ZZ(2**e_phi), sigma_K)

	## Compute [1/q][phi]_*(sigma)(B_push) 
	## = [2**(a-e_phi)/q][sigma]_*(phi)*sigma*phi_dual(B_push_plain)
	# Replaced because it was too slow
	push_sigma_P = a*push_phi(sigma_P)+b*push_phi(sigma_Q)
	push_sigma_Q = c*push_phi(sigma_P)+d*push_phi(sigma_Q)
	# push_phi_P, push_phi_Q = evaluate_isogeny_x_only(push_phi, sigma_P, sigma_Q)
	# push_sigma_P = a*push_phi_P+b*push_phi_Q
	# push_sigma_Q = c*push_phi_P+d*push_phi_Q

	
	return (push_sigma_P, push_sigma_Q)

def update_sig(E_pk, B_Epk_2f_plain, sigma, N, rr, Es, Bs, Kds, e_phi=params.f1):
	r"""
	VUF signature update main function that pushes forward the signature 
	2D-representation through a random 2-isogeny path preomputed with 
	from rerandom_keys.H_path.

	Input:
	- E_pk, public key curve
	- B_Epk_2f_plain, a deterministic basis of E_pk[2^f]
	- sigma = (B_sigma, E_sigma): reformated VUF signature, with
	E_sigma the codomain of sigma: E_pk -> E_sigma of degree N(2^a-N),
	B_sigma the image of (P_pk, Q_pk) via sigma
	- N: prime integer such that deg(sigma) = N(2^a-N).
	- rr: list of random integers between 0 and 2^e_phi-1 determining
	the random 2-isogeny path.
	- Es, Bs, Kds: results from a call to rerandom_keys.H_path (see rerandom_keys.py) 
	using input rr and dual = True.
	- e_phi: integer (default f1, in practice f-f1).

	Output:
	- (B_sigma, E_sigma): updated signature in VUF format.
	"""

	E_start = E_pk
	B_sigma, E_sigma = sigma
	B_start_plain = B_Epk_2f_plain
	cof = 2**(params.f - params.f1)
	B_start = (cof*B_start_plain[0], cof*B_start_plain[1]) 
	B_sigma = (cof*B_sigma[0], cof*B_sigma[1])

	inv_N = inverse_mod(N,2**params.f1)
	B_sigma = (inv_N*B_sigma[0], inv_N*B_sigma[1])

	for i in range(len(rr)):
		# print(f"{i} out of {len(rr)}")
		E_push = Es[i]
		## Compute phi_dual(B_push) with B_push a basis generating E_push[2**f]
		# P_push, Q_push = ec.TorsionBasis(E_push)
		B_push_plain = Bs[i]
		B_push = (cof*B_push_plain[0], cof*B_push_plain[1])
		r = rr[i]
		K_dual = Kds[i]
		

		B_push_sigma = push_forward(E_start,
			E_sigma,B_start_plain,B_start,B_sigma,N,E_push,B_push_plain,r,K_dual,e_phi)
		
		B_push_sigma = (inv_N*B_push_sigma[0], inv_N*B_push_sigma[1])
		
		E_start = E_push
		E_sigma = B_push_sigma[0].curve()
		B_start_plain = B_push_plain
		B_start = B_push
		B_sigma = B_push_sigma

	# Get the codomain Ex
	Ker = ((B_start[0], B_sigma[0]), (B_start[1], B_sigma[1]))
	Phi = hd.Dim2Iso(Ker, params.f1)
	Ex = Phi.codomain()

	return B_start_plain, (B_sigma, Ex)





	

