"""Numerical cross-checks for the three-view derivation, using NumPy.

These checks supplement the written derivations, not an independent review.
Run: python scripts/check-euler-ns-models.py
"""
import numpy as np

rng = np.random.default_rng(20261109)
R, cv = 1.0, 2.5
cp, gamma = cv + R, 1 + R / cv


def unpack(U):
    rho, momentum, energy = U
    u = momentum / rho
    T = (energy / rho - u * u / 2) / cv
    p = rho * R * T
    s = cv * np.log(T) - R * np.log(rho)
    return rho, u, T, p, s


def flux(U):
    rho, u, T, p, s = unpack(U)
    return np.array([U[1], U[1] * u + p, (U[2] + p) * u])


def eta(U):
    rho, u, T, p, s = unpack(U)
    return -rho * s


def entropy_flux(U):
    rho, u, T, p, s = unpack(U)
    return -rho * s * u


def entropy_variables(U):
    rho, u, T, p, s = unpack(U)
    return np.array([cp - s - u * u / (2 * T), u / T, -1 / T])


def jacobian(fun, U):
    step = 1e-25
    return np.array([
        np.imag(fun(U.astype(complex) + 1j * step * e)) / step
        for e in np.eye(3)
    ]).T


for _ in range(60):
    rho, T = np.exp(rng.uniform(-1, 1, 2))
    u = rng.uniform(-3, 3)
    c = np.sqrt(gamma * R * T)
    p = rho * R * T
    H = cp * T + u * u / 2
    U = np.array([rho, rho * u, rho * (cv * T + u * u / 2)])
    AW = np.array([[u, rho, 0], [0, u, 1/rho], [0, rho*c*c, u]])
    AU = np.array([
        [0, 1, 0],
        [(gamma-3)*u*u/2, (3-gamma)*u, gamma-1],
        [u*((gamma-1)*u*u/2-H), H-(gamma-1)*u*u, gamma*u],
    ])
    B = np.array([[1,0,0],[u,rho,0],[u*u/2,rho*u,1/(gamma-1)]])
    RW = np.array([[1,1,1],[-c/rho,0,c/rho],[c*c,0,c*c]])
    LW = np.array([[0,-rho/(2*c),1/(2*c*c)],[1,0,-1/(c*c)],[0,rho/(2*c),1/(2*c*c)]])
    RU = np.array([[1,1,1],[u-c,u,u+c],[H-u*c,u*u/2,H+u*c]])
    lam = np.diag([u-c,u,u+c])
    np.testing.assert_allclose(AU, jacobian(flux,U), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(AU@B, B@AW, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(AW@RW, RW@lam, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(LW@RW, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(AU@RU, RU@lam, rtol=1e-12, atol=1e-12)
    for ell in rng.uniform(-4,4,3):
        np.testing.assert_allclose(np.linalg.det(AU-ell*np.eye(3)), (u-ell)*((u-ell)**2-c*c), rtol=1e-10, atol=1e-10)
    AL=np.array([[0,-1,0],[-rho*rho*c*c,0,p/cv],[0,0,0]])
    np.testing.assert_allclose(np.sort(np.linalg.eigvals(AL)),[-rho*c,0,rho*c],atol=1e-11)
    np.testing.assert_allclose(u+np.sort(np.linalg.eigvals(AL))/rho,[u-c,u,u+c],atol=1e-11)

    z=entropy_variables(U)
    hess=jacobian(entropy_variables,U)
    np.testing.assert_allclose(jacobian(eta,U),z,rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(jacobian(entropy_flux,U),z@AU,rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(hess,hess.T,atol=1e-11)
    assert np.min(np.linalg.eigvalsh(hess))>0
    np.testing.assert_allclose(hess@AU,(hess@AU).T,atol=2e-11)
    C=np.array([[1,0,0],[u,rho,0],[cv*T+u*u/2,rho*u,rho*cv]])
    np.testing.assert_allclose(C.T@hess@C,np.diag([R/rho,rho/T,rho*cv/T**2]),rtol=1e-11,atol=1e-11)
print('PASS: 60 states, flux derivatives, similarities, characteristic polynomial and eigenvectors, mass-coordinate speeds')
print('PASS: entropy variables, flux compatibility, positive Hessian and symmetrization')

for _ in range(40):
    rho,T=np.exp(rng.uniform(-1,1,2)); u=rng.uniform(-2,2)
    mu,zeta,k=np.exp(rng.uniform(-1,1,3))
    av=4*mu/3+zeta
    rhox,ux,Tx=rng.normal(size=3)
    q=-k*Tx; tau=av*ux
    zx=np.array([-cv*Tx/T+R*rhox/rho-u*ux/T+u*u*Tx/(2*T*T),ux/T-u*Tx/(T*T),Tx/(T*T)])
    Fv=np.array([0,tau,u*tau-q])
    np.testing.assert_allclose(zx@Fv,av*ux*ux/T+k*Tx*Tx/(T*T))
    U=np.array([rho,rho*u,rho*(cv*T+u*u/2)])
    np.testing.assert_allclose(entropy_variables(U)@Fv,q/T,atol=1e-12)
    nu,chi=av/rho,k/(rho*cv)
    c2=gamma*R*T
    A=np.array([[u,rho,0],[R*T/rho,u,R],[0,(gamma-1)*T,u]])
    B=np.diag([0,nu,chi])
    weight=np.diag([R/rho,rho/T,rho*cv/T**2])
    np.testing.assert_allclose(weight@A,(weight@A).T,atol=1e-12)
    for wave in [0.02,0.3,1.7]:
        M=-1j*wave*A-wave*wave*B
        omega=np.linalg.eigvals(M)
        assert max(omega.real)<1e-11
        roots=np.roots([1,(nu+chi)*wave**2,c2*wave**2+nu*chi*wave**4,c2*chi*wave**4/gamma])
        for root in roots:
            assert min(abs(root-(omega+1j*wave*u)))<1e-9
        A5=np.eye(5)*u
        block=[0,1,4]
        A5[np.ix_(block,block)]=A
        B5=np.diag([0,nu,mu/rho,mu/rho,chi])
        np.testing.assert_allclose(np.sort(np.linalg.eigvals(A5)),[u-np.sqrt(c2),u,u,u,u+np.sqrt(c2)],atol=1e-12)
        roots5=np.linalg.eigvals(-1j*wave*A5-wave*wave*B5)
        shear=-1j*wave*u-mu/rho*wave*wave
        assert np.count_nonzero(abs(roots5-shear)<1e-10)>=2
        for root in omega:
            assert min(abs(root-roots5))<1e-9
    wave=1e-4
    roots=np.linalg.eigvals(-1j*wave*A-wave*wave*B)
    DT=chi/gamma
    for speed,damping in [(u,DT),(u-np.sqrt(c2),(nu+(gamma-1)*DT)/2),(u+np.sqrt(c2),(nu+(gamma-1)*DT)/2)]:
        target=-1j*wave*speed-damping*wave**2
        assert min(abs(roots-target))<2e-10
print('PASS: viscous entropy flux/production and 40 constant-state Fourier spectra, including long-wave limits and 3D shear modes')

for _ in range(30):
    rotation,_=np.linalg.qr(rng.normal(size=(3,3)))
    if np.linalg.det(rotation)<0:
        rotation[:,0]*=-1
    F=rotation@np.diag(np.exp(rng.uniform(-.5,.5,3)))
    J=np.linalg.det(F); inv=np.linalg.inv(F)
    gradV=rng.normal(size=(3,3)); L=gradV@inv
    sigma=rng.normal(size=(3,3)); sigma=(sigma+sigma.T)/2
    P=J*sigma@inv.T
    V,q,N=rng.normal(size=(3,3))
    Q=J*inv@q; area=J*inv.T@N
    np.testing.assert_allclose(np.sum(P*gradV),J*np.sum(sigma*L),atol=1e-12)
    np.testing.assert_allclose(Q@N,q@area,atol=1e-12)
    np.testing.assert_allclose((P@N)@V,(sigma@area)@V,atol=1e-12)
    gradT=rng.normal(size=3); k=1.7
    np.testing.assert_allclose(-J*k*inv@inv.T@(F.T@gradT),J*inv@(-k*gradT),atol=1e-12)
print('PASS: reference stress power, Nanson traction and heat-flux transformations')

g=1.4; rhoL=pL=1.; uL=2*np.sqrt(g)
rhoR,pR=8/3,4.5; uR=uL/rhoR
UL=np.array([rhoL,rhoL*uL,pL/(g-1)+rhoL*uL*uL/2])
UR=np.array([rhoR,rhoR*uR,pR/(g-1)+rhoR*uR*uR/2])
np.testing.assert_allclose(flux(UL),flux(UR),atol=1e-12)
assert np.log(pR/pL)-g*np.log(rhoR/rhoL)>0
print('PASS: Mach-2 stationary shock mass/momentum/energy jumps and entropy increase')
