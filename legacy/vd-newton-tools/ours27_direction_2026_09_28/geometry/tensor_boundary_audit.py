"""Tensor-scale improper posterior limits: proper, but cutoff-dependent laws.

Scratch audit. No tensor candidate has been adopted. Managed context is
VD-docs@b4bfc2040031ef010416334d0c5f4aa2ad05cbd8.

SETUP
-----
For R positive definite of dimension d, take q(R)=|I+R|^{-(d+1)/2} with
Lebesgue measure on symmetric matrices. With N readings, and S=sum A_i A_i',
the posterior after Beta=R(I+R)^(-1) is

  |I-Beta|^alpha exp(tr(Beta*S)/2) dBeta,
  alpha=(N-d-1)/2, 0<Beta<I.

Indeed dR=|I-Beta|^{-(d+1)}dBeta. Put X=I-Beta; the real symmetric eigenvalue
Jacobian gives eigen-deficit density prod(x_i^alpha)*prod_{i<j}|x_i-x_j|,
times Haar orientation. The bounded, positive data tilt can be restored to
any weak limit by ordinary normalization. All conclusions below therefore
extend from S=0 to each fixed finite S.

A cluster of k deficits of order delta contributes the power
  delta^[k(alpha+1)+k(k-1)/2] = delta^[k(N-d+k)/2].
This includes the Vandermonde repulsion. Counting only individual divergent
factors incorrectly predicts how many eigenvalues reach the boundary.

d=2, N=1: ONE LOGARITHMIC ESCAPE
--------------------------------
Here alpha=-1. Under the spectral ceiling R<L I, every deficit exceeds
delta=1/(L+1). For ordered x<y, the density is (y-x)/(x*y). Integrating x
gives log(1/delta)+O(1) for every fixed y>0. The normalized boundary law is
  X = y*v*v', y~Uniform(0,1), v uniform on the unoriented unit circle.
Thus Beta has one eigenvalue 1 and one uniform eigenvalue in (0,1).
Any fixed positive ellipsoidal ceiling R<L C changes log L only by a bounded
orientation term, so this particular logarithmic limit is unchanged.

d=3, N=1: ORDINARY SPECTRAL CEILING HAS ONLY ONE ESCAPING EIGENVALUE
----------------------------------------------------------------------------
Now alpha=-3/2. Write ordered deficits 0<x1<x<y<1. Before integration their
kernel is
  x1^(-3/2)*(xy)^(-3/2)*(x-x1)*(y-x1)*(y-x).
Under x1>delta, integration gives leading
  2*delta^(-1/2)*(xy)^(-1/2)*(y-x).
The remaining normalized ordered-deficit law is
  p_spectral(x,y)=(3/2)*(xy)^(-1/2)*(y-x), 0<x<y<1.
It is proper. Equivalently y~Beta(2,1), x/y~Beta(1/2,2), independently.
Consequently E tr(X)=4/5 and E tr(Beta)=11/5=2.2. X has rank TWO, so only
ONE eigenvalue of Beta equals 1. Dividing the finite-cutoff x1 integral by
2*delta^(-1/2) bounds it above by (xy)^(-1/2)*(y-x), an integrable function;
dominated convergence supplies the full normalized boundary limit.

ORTHOGONALLY INVARIANT DETERMINANT CEILING GIVES A DIFFERENT LIMIT
----------------------------------------------------------------
Instead truncate |I+R|<L, or x1*x*y>1/L. This also exhausts the entire cone
and is invariant under every orthogonal recoding. Now x1>1/(L*x*y), so its
leading integrated factor is 2*sqrt(L*x*y). The remaining kernel becomes
  y-x.
After normalization,
  p_det(x,y)=6*(y-x), 0<x<y<1.
Equivalently y~Beta(3,1), x/y~Beta(1,2), independently. E tr(X)=1 and
E tr(Beta)=2. Dividing by 2*sqrt(L), the exact integrated density is bounded
above by y-x, again proving the limit by dominated convergence. Therefore
even orthogonal invariance and exhaustion do not select a unique boundary
probability law for this improper hyperposterior.

ANISOTROPIC HOMOTHETIC CEILINGS ALSO LEAVE A TRACE
------------------------------------------------
For d=3,N=1 and R<L C, the smallest deficit along the escaping direction w
has leading lower bound (w'C^(-1)w)/L. Its limit therefore weights the Haar
orientation by (w'C^(-1)w)^(-1/2). For C=diag(4,1,1),
  E[w_1^2]=2/3-sqrt(3)/(2*pi)=0.391002..., rather than 1/3.
This contrasts with the logarithmic d=2 case, where fixed anisotropy changes
only a lower-order additive constant.

SINGULAR WISHART / MATRIX-BETA EXTENSION
---------------------------------------
There is another proper construction: A~Wishart_d(d+1,I), B~Wishart_d(N,I),
and Beta=(A+B)^(-1/2) A (A+B)^(-1/2). For N>=d its ordinary matrix-Beta
density is exactly the displayed kernel without its data tilt. For N<d,
I-Beta has rank N, so Beta has d-N eigenvalues exactly 1.

When N=1, I-Beta=x*v*v' with uniform orientation and x~Beta(d/2,1). This
follows from the squared row norm of a d-frame in R^(d+2), equivalently the
diagonal leverage of a uniform rank-d projection. Thus this extension agrees
with the d=2,N=1 spectral limit, but for d=3,N=1 it has E tr(Beta)=12/5=2.4
and rank(I-Beta)=1. It is different from BOTH cutoff limits above. It is a
well-defined proper continuation, not the normalized limit of the stated
spectral truncation. A bounded exponential data tilt preserves propriety.

Primary-source check: Diaz-Garcia and Gutierrez-Jaimez, arXiv:0904.2138,
Definition 2.1 and equations (1), (3), https://arxiv.org/pdf/0904.2138.
For rank q=r=1 and s=d+1, the rank-one density contributes x^(-d/2) and
the Hausdorff eigenvalue element contributes x^(d-1), giving x^(d/2-1),
exactly the Beta(d/2,1) scalar density above. The cutoff comparison and
sequential-limit obstruction here are our direct calculations, not claims
attributed to that paper.

PROPER CURRENT CANDIDATES DO NOT RESTORE ORDINARY SEQUENTIAL UPDATING
-------------------------------------------------------------------
Every boundary law above lives on the compact matrix interval 0<=Beta<=I.
For each current standardized candidate, v_i|Beta,A_i is the proper Gaussian
N(Beta*A_i,Beta), including its ordinary or degenerate Gaussian boundary
interpretation. Thus escaping R does not by itself prevent proper current
candidate laws.

However, the next reading's Gaussian marginal likelihood contains
  |I-Beta|^(1/2)*exp(tr(Beta*A_new*A_new')/2)
(times a finite data-only factor). It is identically zero on these boundary
supports. Multiplying the normalized N-reading boundary law by this factor
cannot produce an (N+1)-reading posterior. For the singular-Wishart extension
the support also changes rank with N; ordinary likelihood updating cannot
create mass outside the old support. Finite cutoffs can be updated before
taking their limits, but normalization/limit and updating do not commute.
The new-reading predictive mass escapes in the directions of infinite R.
This is a genuine coherence issue if ordinary sequential Bayes updating is
required, not a claim that a batch-defined operational rule is meaningless.
"""
from pathlib import Path
from fractions import Fraction
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss


def finite_cutoff_mean_trace_x(kind,L,order=200):
    """Integrate the smallest deficit analytically, the other two by quadrature."""
    z,w=leggauss(order);q=(z+1)/2;w=w/2
    delta=1/(L+1) if kind=='spectral' else 1/L
    ymin=delta if kind=='spectral' else delta**(1/3)
    # Square-root coordinates remove the limiting spectral residual singularity.
    sy=np.sqrt(ymin)+(1-np.sqrt(ymin))*q
    y=sy*sy;wy=w*(1-np.sqrt(ymin))*2*sy
    xmin=np.full_like(y,delta) if kind=='spectral' else np.sqrt(delta/y)
    sx=np.sqrt(xmin)[:,None]+(np.sqrt(y)-np.sqrt(xmin))[:,None]*q[None,:]
    x=sx*sx
    wx=w[None,:]*(np.sqrt(y)-np.sqrt(xmin))[:,None]*2*sx
    yy=y[:,None]
    lo=np.full_like(x,delta) if kind=='spectral' else delta/(x*yy)
    def anti0(t):
        return -2*x*yy/np.sqrt(t)-2*(x+yy)*np.sqrt(t)+(2/3)*t**1.5
    def anti1(t):
        return 2*x*yy*np.sqrt(t)-(2/3)*(x+yy)*t**1.5+(2/5)*t**2.5
    i0=np.maximum(0.,anti0(x)-anti0(lo))
    i1=np.maximum(0.,anti1(x)-anti1(lo))
    rest=(x*yy)**(-1.5)*(yy-x)
    weights=wy[:,None]*wx*rest
    normalizer=float(np.sum(weights*i0))
    trace=float(np.sum(weights*(i1+(x+yy)*i0))/normalizer)
    return {'cutoff':kind,'L':L,'normalizer_without_Haar_constant':normalizer,
            'mean_trace_X':trace,'mean_trace_Beta':3-trace}


def main():
    rows=[]
    for L in (1e3,1e5,1e7,1e9):
        for kind in ('spectral','determinant'):
            rows.append(finite_cutoff_mean_trace_x(kind,L))
    assert abs(rows[-2]['mean_trace_X']-.8)<.002
    assert abs(rows[-1]['mean_trace_X']-1)<.004
    out={
        'source_commit':'b4bfc2040031ef010416334d0c5f4aa2ad05cbd8',
        'd3_N1_limits':[
            {'construction':'spectral ceiling','rank_I_minus_Beta':2,
             'E_trace_X':float(Fraction(4,5)),'E_trace_Beta':float(Fraction(11,5))},
            {'construction':'determinant ceiling','rank_I_minus_Beta':2,
             'E_trace_X':1.,'E_trace_Beta':2.},
            {'construction':'singular Wishart ratio','rank_I_minus_Beta':1,
             'E_trace_X':float(Fraction(3,5)),'E_trace_Beta':float(Fraction(12,5))}],
        'd2_N1_spectral_and_singular_Wishart_agree':True,
        'anisotropic_spectral_C_diag_4_1_1_E_escape_w1_squared':2/3-math.sqrt(3)/(2*math.pi),
        'isotropic_E_escape_w1_squared':1/3,
        'finite_cutoff_quadrature':rows,
        'quadrature_order_200_to_400_mean_trace_differences_at_L_1e9':{
            kind:abs(finite_cutoff_mean_trace_x(kind,1e9,200)['mean_trace_X']-
                     finite_cutoff_mean_trace_x(kind,1e9,400)['mean_trace_X'])
            for kind in ('spectral','determinant')},
        'primary_source_singular_beta_rank_one_check':'https://arxiv.org/pdf/0904.2138, Definition 2.1 and equations (1), (3)',
        'status':'Proper current-candidate boundary laws exist, but are cutoff-dependent and do not support ordinary post-limit sequential Bayes updating.'}
    Path(__file__).with_name('tensor_boundary_audit_results.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
