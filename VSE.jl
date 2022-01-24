
const me = 0.511 # MeV
const ec = 0.302811
const αem = ec^2/(4pi)

#################################################
# VECTOR SELF-ENERGY IN DILUTE NON-RELATIVISTIC GAS

# Plasma dispersion relation Z(xi) in a form valid for real arguments
function Zi(x)
	sqrt(pi)*exp(-x^2)
end

function Zr(x)
	-2*dawson(x)
end

# Real part of longitudinal photon polarization tensor in dilute, non-relativistic gas.
function ΠLr_dnr(q0,q,m,Z,n,T)
	s = sqrt(T/m)
	xi = q0/(sqrt(2)*q*s)
	qp = (q0^2 - q^2) / (2*q0*m)
	(
	 ((ec^2 * Z^2 * n) / (sqrt(2)*q*s)) *
	 (Zr(xi*(1 + qp)) - Zr(xi*(1 - qp)))
	 )
end

# Real part of transverse photon polarization tensor in dilute, non-relativistic gas.
function ΠTr_dnr(q0,q,m,Z,n,T)
	wp2 = ec^2 * Z^2 * n / m
	s = sqrt(T/m)
	xi = q0/(sqrt(2)*q*s)
	Q2 = q0^2 - q^2
	qp = Q2 / (2*q0*m)
	(
	 wp2*(1 - (s^2 - Q2/(4*m^2))*(m/q0)*xi*
		  (Zr(xi*(1 + qp)) - Zr(xi*(1 - qp))))
	 )
end

# Imaginary part of longitudinal photon polarization tensor in dilute, non-relativistic gas
function ΠLi_dnr(q0,q,m,Z,n,T)
	s = sqrt(T/m)
	xi = q0/(sqrt(2)*q*s)
	qp = (q0^2 - q^2) / (2*q0*m)
	(
	 ((ec^2 * Z^2 * n) / (sqrt(2)*q*s)) *
	 (Zi(xi*(1 + qp)) - Zi(xi*(1 - qp)))
	 )
end

# Imaginary part of transverse photon polarization tensor in dilute, non-relativistic gas
function ΠTi_dnr(q0,q,m,Z,n,T)
	wp2 = ec^2 * Z^2 * n / m
	s = sqrt(T/m)
	xi = q0/(sqrt(2)*q*s)
	Q2 = q0^2 - q^2
	qp = Q2 / (2*q0*m)
	(
	 -wp2*((s^2 - Q2/(4*m^2))*(m/q0)*xi*
		  (Zi(xi*(1 + qp)) - Zi(xi*(1 - qp))))
	 )
end

#################################################
# VECTOR SELF-ENERGY IN DEGENERATE FERMION GAS 

# Imaginary part of longitudinal photon polarization tensor for degenerate electron gas.
function ΠLi_e_degen(q0,q,mue;mf=me)
	Q2 = q0^2 - q^2
	pre = -(ec^2/(4pi)) * (Q2/q^3)
	ep = max(mue,mf + q0,0.5*(q0 + q*sqrt(1 - 4*mf^2/Q2)))
	if ep > mue + q0
		0.
	else
		pre * (-((2*ep^3)/3) + (2*mue^3)/3 + ep^2*q0 + mue^2*q0 - q0^3/3 
		 - (ep*Q2)/2 + (mue*Q2)/2 + (q0*Q2)/2)
	end
end

# Imaginary part of transverse photon polarization tensor for degenerate electron gas.
function ΠTi_e_degen(q0,q,mue;mf=me)
	Q2 = q0^2 - q^2
	pre = -(ec^2/(4pi)) * (1/q)
	ep = max(mue,mf + q0,0.5*(q0 + q*sqrt(1 - 4*mf^2/Q2)))
	if ep > mue + q0
		0.
	else
		(pre/12)*(-12*mf^2*(-ep + mue + q0) + 4*(-ep^3 + (mue + q0)^3) - 
 (4*q0^2*(-ep^3 + (mue + q0)^3))/q^2 - 6*(-ep + mue + q0)*Q2 + 
 (6*q0*(-ep^2 + (mue + q0)^2)*Q2)/q^2 - (3*(-ep + mue + q0)*Q2^2)/q^2)
	end
end

# Real part of longitudinal photon polarization tensor for degenerate electron gas.
# For internal integral definitions, see below.
function ΠLr_e_degen(q0,q,mue;mf=me)
	Q2 = q0^2 - q^2

	pf = sqrt(mue^2 - mf^2)

	ik = q*(-(mue*pf) + mf^2*log((mue + pf)/mf))

	i0m = integ0m(q0,q,mue,pf,mf,Q2)
	i1m = integ1m(q0,q,mue,pf,mf,Q2)
	i2m = integ2m(q0,q,mue,pf,mf,Q2)

	i0p = integ0p(q0,q,mue,pf,mf,Q2)
	i1p = integ1p(q0,q,mue,pf,mf,Q2)
	i2p = integ2p(q0,q,mue,pf,mf,Q2)

	int = (ik - (Q2/4)*i0m - i2m + q0*i1m
		   + i2p + (Q2/4)*i0p + q0*i1p)
	(Q2/q^3)*(ec^2/(2*pi^2))*int
end

# Real part of transverse photon polarization tensor for degenerate electron gas.
# For internal integral definitions, see below.
function ΠTr_e_degen(q0,q,mue;mf=me)
	Q2 = q0^2 - q^2

	pf = sqrt(mue^2 - mf^2)

	ik = -(1 + Q2/(2*q^2))*q*(-(mue*pf) + mf^2*log((mue + pf)/mf))

	i0m = integ0m(q0,q,mue,pf,mf,Q2)
	i1m = integ1m(q0,q,mue,pf,mf,Q2)
	i2m = integ2m(q0,q,mue,pf,mf,Q2)

	i0p = integ0p(q0,q,mue,pf,mf,Q2)
	i1p = integ1p(q0,q,mue,pf,mf,Q2)
	i2p = integ2p(q0,q,mue,pf,mf,Q2)

	int = (ik + (Q2/4 + Q2^2/(8*q^2) + mf^2/2)*(i0m - i0p) 
		   + (Q2/(2*q^2))*(i2m-i2p) - q0*Q2*(i1m + i1p)/(2*q^2))
	(ec^2/(2*pi^2))*int/q
end

# Helper functions and internal integrals for real parts of polarization tensors in degenerate electron gas.

function logabs(x)
	log(abs(x))
end

function arcTanh(x)
	0.5*log(abs((1 + x)/(1 - x)))
end

function integ0m(q0,q,mue,pf,mf,Q2)
	QQ2 = q^2 + q0^2
	C1 = sqrt(Q2*(Q2 - 4*mf^2))
	((2*q*logabs((mue + pf)/mf) + 
	 2*mue*logabs((-2*pf*q - 2*mue*q0 + Q2)/(2*pf*q - 2*mue*q0 + Q2)) + 
	 ((4*mf^2*q*q0 - 2*q*q0*Q2 + C1*QQ2)*(logabs(mf) + logabs(-(C1*q) - 2*mue*Q2 + q0*Q2)))/
	 sqrt(Q2*(-2*q0*(C1*q + 2*mf^2*q0) + Q2*QQ2)) + 
	 ((2*q*q0*(-2*mf^2 + Q2) + C1*QQ2)*(-logabs(mf) - logabs(C1*q - 2*mue*Q2 + q0*Q2)))/
	 sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)) - 
	 ((4*mf^2*q*q0 - 2*q*q0*Q2 + C1*QQ2)*logabs(C1*mue*q + 2*mf^2*Q2 - mue*q0*Q2 + 
											 pf*sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))))/
	 sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)) + 
	 ((-4*mf^2*q*q0 + 2*q*q0*Q2 + C1*QQ2)*logabs(C1*mue*q - 2*mf^2*Q2 + mue*q0*Q2 - 
											  pf*sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))))/
	 sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))/2)
end

function integ0p(q0,q,mue,pf,mf,Q2)
	QQ2 = q^2 + q0^2
	C1 = sqrt(Q2*(Q2 - 4*mf^2))

	((-2*q*logabs((mue + pf)/mf) + 
  2*mue*logabs((2*pf*q + 2*mue*q0 + Q2)/(-2*pf*q + 2*mue*q0 + Q2)) + 
  ((4*mf^2*q*q0 - 2*q*q0*Q2 - C1*QQ2)*(logabs(mf) + logabs(-(C1*q) - (2*mue + q0)*Q2)))/
   sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)) - 
  ((4*mf^2*q*q0 - 2*q*q0*Q2 + C1*QQ2)*(-logabs(mf) - logabs(C1*q - (2*mue + q0)*Q2)))/
   sqrt(Q2*(-2*q0*(C1*q + 2*mf^2*q0) + Q2*QQ2)) - 
  ((4*mf^2*q*q0 - 2*q*q0*Q2 + C1*QQ2)*logabs(C1*mue*q - 2*mf^2*Q2 - mue*q0*Q2 - 
      pf*sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))))/
   sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)) + 
  ((-4*mf^2*q*q0 + 2*q*q0*Q2 + C1*QQ2)*logabs(C1*mue*q + 2*mf^2*Q2 + mue*q0*Q2 + 
      pf*sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))))/
  sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))/2)
end

function integ1m(q0,q,mue,pf,mf,Q2)
	QQ2 = q^2 + q0^2
	C1 = sqrt(Q2*(Q2 - 4*mf^2))
	((pf*q + (q*q0*(-2*mf^2 + Q2)*arcTanh(mue/pf))/Q2 - 
		((8*q0*Q2*(-2*mf^2 + Q2)*(4*mf^2*q^2 + Q2^2) + (4*C1*q - 4*q0*Q2)*
		  (4*mf^2*(q^4 + q^2*q0^2 - 2*q0^4) + (q^2 + 3*q0^2)*Q2^2))*
		 sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))*
		 arcTanh((-(C1*mue*q) - 2*mf^2*Q2 + mue*q0*Q2)/(sqrt(-mf^2 + mue^2)*
														sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))))/
		(C1*Q2*(-64*mf^2*Q2^2 + (4*C1*q - 4*q0*Q2)^2)) - 
		((16*mf^4*q0*(q^3 - q*q0^2)^2 + Q2^3*(C1*(q^3 + 3*q*q0^2) + q0*(3*q^2 + q0^2)*Q2) - 
		  4*mf^2*Q2^2*(C1*(q^3 + 2*q*q0^2) + q0*(4*q^2 + q0^2)*Q2))*
		 arcTanh((C1*mue*q - 2*mf^2*Q2 + mue*q0*Q2)/(sqrt(-mf^2 + mue^2)*
													 sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))))/
		(4*C1*Q2^2*sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))) + 
		mue^2*logabs((-2*pf*q - 2*mue*q0 + Q2)/(2*pf*q - 2*mue*q0 + Q2)))/2)
end

function integ1p(q0,q,mue,pf,mf,Q2)
	QQ2 = q^2 + q0^2
	C1 = sqrt(Q2*(Q2 - 4*mf^2))
	((-(pf*q) + (q*q0*(-2*mf^2 + Q2)*arcTanh(mue/pf))/Q2 + 
  ((8*q0*Q2*(-2*mf^2 + Q2)*(4*mf^2*q^2 + Q2^2) - (-4*C1*q + 4*q0*Q2)*
      (4*mf^2*(q^4 + q^2*q0^2 - 2*q0^4) + (q^2 + 3*q0^2)*Q2^2))*
    sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))*
    arcTanh((C1*mue*q - (2*mf^2 + mue*q0)*Q2)/
      (pf*sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))))/
   (C1*Q2*(-64*mf^2*Q2^2 + (4*C1*q - 4*q0*Q2)^2)) + 
  ((16*mf^4*q0*(q^3 - q*q0^2)^2 + Q2^3*(C1*(q^3 + 3*q*q0^2) + q0*(3*q^2 + q0^2)*Q2) - 
     4*mf^2*Q2^2*(C1*(q^3 + 2*q*q0^2) + q0*(4*q^2 + q0^2)*Q2))*
    arcTanh((-2*mf^2*Q2 - mue*(C1*q + q0*Q2))/
      (pf*sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))))/
   (4*C1*Q2^2*sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))) + 
  mue^2*logabs((2*pf*q + 2*mue*q0 + Q2)/(-2*pf*q + 2*mue*q0 + Q2)))/2)
end

function integ2m(q0,q,mue,pf,mf,Q2)
	QQ2 = q^2 + q0^2
	C1 = sqrt(Q2*(Q2 - 4*mf^2))
	((2*mue*pf*q + (4*pf*q*q0*(-2*mf^2 + Q2))/Q2 + 
  2*mf^2*q*arcTanh(mue/pf) + 
  (q*(4*mf^2*(q^4 + q^2*q0^2 - 2*q0^4) + (q^2 + 3*q0^2)*Q2^2)*
    arcTanh(mue/pf))/Q2^2 - 
  (8*sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))*
    (Q2^2*(4*mf^2*q^2 + Q2^2)*(4*mf^2*(q^4 + q^2*q0^2 - 2*q0^4) + 
       (q^2 + 3*q0^2)*Q2^2) + q0*(4*C1*q - 4*q0*Q2)*(4*mf^4*(q^3 - q*q0^2)^2 + 
       mf^2*(5*q^4 - 2*q^2*q0^2 - 3*q0^4)*Q2^2 + Q2^4*QQ2))*
    arcTanh((-(C1*mue*q) - 2*mf^2*Q2 + mue*q0*Q2)/(sqrt(-mf^2 + mue^2)*
       sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))))/
   (C1*Q2^3*(-64*mf^2*Q2^2 + (4*C1*q - 4*q0*Q2)^2)) + 
  (sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))*
    (Q2^2*(4*mf^2*q^2 + Q2^2)*(4*mf^2*(q^4 + q^2*q0^2 - 2*q0^4) + 
       (q^2 + 3*q0^2)*Q2^2) - 4*q0*(C1*q + q0*Q2)*(4*mf^4*(q^3 - q*q0^2)^2 + 
       mf^2*(5*q^4 - 2*q^2*q0^2 - 3*q0^4)*Q2^2 + Q2^4*QQ2))*
    arcTanh((C1*mue*q - 2*mf^2*Q2 + mue*q0*Q2)/(sqrt(-mf^2 + mue^2)*
       sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))))/
   (2*C1*Q2^3*(-4*mf^2*Q2^2 + (C1*q + q0*Q2)^2)) + 
   4*mue^3*logabs((-2*pf*q - 2*mue*q0 + Q2)/(2*pf*q - 2*mue*q0 + Q2)))/12)
end

function integ2p(q0,q,mue,pf,mf,Q2)
	QQ2 = q^2 + q0^2
	C1 = sqrt(Q2*(Q2 - 4*mf^2))
	((-2*mue*pf*q + (4*pf*q*q0*(-2*mf^2 + Q2))/Q2 - 2*mf^2*q*arcTanh(mue/pf) - 
  (q*(4*mf^2*(q^4 + q^2*q0^2 - 2*q0^4) + (q^2 + 3*q0^2)*Q2^2)*arcTanh(mue/pf))/Q2^2 - 
  (sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))*
    (8*Q2^2*(4*mf^2*q^2 + Q2^2)*(4*mf^2*(q^4 + q^2*q0^2 - 2*q0^4) + 
       (q^2 + 3*q0^2)*Q2^2) - 8*q0*(-4*C1*q + 4*q0*Q2)*(4*mf^4*(q^3 - q*q0^2)^2 + 
       mf^2*(5*q^4 - 2*q^2*q0^2 - 3*q0^4)*Q2^2 + Q2^4*QQ2))*
    arcTanh((C1*mue*q - (2*mf^2 + mue*q0)*Q2)/
      (pf*sqrt(Q2*(-2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))))/
   (C1*Q2^3*(-64*mf^2*Q2^2 + (4*C1*q - 4*q0*Q2)^2)) + 
  (sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2))*
    (8*Q2^2*(4*mf^2*q^2 + Q2^2)*(4*mf^2*(q^4 + q^2*q0^2 - 2*q0^4) + 
       (q^2 + 3*q0^2)*Q2^2) - 32*q0*(C1*q + q0*Q2)*(4*mf^4*(q^3 - q*q0^2)^2 + 
       mf^2*(5*q^4 - 2*q^2*q0^2 - 3*q0^4)*Q2^2 + Q2^4*QQ2))*
    arcTanh((-2*mf^2*Q2 - mue*(C1*q + q0*Q2))/
      (pf*sqrt(Q2*(2*C1*q*q0 - 4*mf^2*q0^2 + Q2*QQ2)))))/
   (16*C1*Q2^3*(-4*mf^2*Q2^2 + (C1*q + q0*Q2)^2)) + 
  4*mue^3*logabs((2*pf*q + 2*mue*q0 + Q2)/(-2*pf*q + 2*mue*q0 + Q2)))/12)
end
