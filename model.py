import numpy as np
from numpy import sum  
from numpy.random import normal as gran
from numpy import pi as π

def dV(x, param):
    ωc = param.ωc
    χ =  param.χ
    dv = np.zeros(len(x))
    c = χ * (2.0 /ωc**3.0)**0.5
    #-----------------------------------------------------------------!
    # Nuclear DOFs ---------------------------------------------------!
    #-----------------------------------------------------------------!
    dv[0]  =  dE(x[0])                           # Molecular part
    dv[0] +=  (ωc**2.0) * (x[1] + c * x[0]) * c  # Cavity

    ωbj, cbj =  param.ωbj, param.cbj
    nb = len(ωbj)
    # Bath
    dv[0] += sum((ωbj**2.0) * (x[2:2+nb] + cbj * x[0]/ωbj**2.0) * (cbj/ωbj**2.0))  
    #-----------------------------------------------------------------!
    # Cavity DOFs ----------------------------------------------------!
    #-----------------------------------------------------------------!
    dv[1]   =  (ωc**2.0) * (x[1] + c * x[0])  
    ωcj, ccj =  param.ωcj, param.ccj
    dv[1] += sum((ωcj**2.0) * (x[2+nb:] + ccj * x[1]/ωcj**2.0) * (ccj/ωcj**2.0)) 
    #-----------------------------------------------------------------!
    # Solvent DOFs ---------------------------------------------------!
    #-----------------------------------------------------------------!
    dv[2:2+nb] = (ωbj**2.0) * (x[2:2+nb] + cbj * x[0]/ωbj**2.0 )
    #-----------------------------------------------------------------!
    # Loss DOFs    ---------------------------------------------------!
    #-----------------------------------------------------------------!
    dv[2+nb:] = (ωcj**2.0) * (x[2+nb:] + ccj * x[1]/ωcj**2.0) 
    return dv

def force(x, param):
    return -1.0 * dV(x, param)

def E(R):
    ωb = 0.004556335 # 1000 cm-1
    Eb = 0.01025175  # 2250 cm-1
    a =  ωb**4/(16*Eb) #/10.0
    b = -ωb**2/2 #/10.0

    return a*(R**4) + b*(R**2) #+ c * R

def dE(R):
    ωb = 0.004556335 # 1000 cm-1
    Eb = 0.01025175  # 2250 cm-1
    a =  ωb**4/(16*Eb) #/10.0
    b = -ωb**2/2 #/10.0
    return 4.0 * a * (R**3) + 2.0 * b * R #+ c * R

 
def init(param):
    m = param.m
    β = param.β
    ndof = param.ndof

    # Cavity
    ωc = param.ωc
    χ = param.χ
    # cavity Loss
    ωcj, ccj = param.ωcj, param.ccj
    # Solvent
    ωbj, cbj = param.ωbj, param.cbj

    # All mode frequencies
    ω = np.zeros((ndof))
    ω[1] = param.ωc              # 1 --> Cavity
    ω[2:param.nbath+2] = param.ωbj # 2: --> Solvent
    ω[param.nbath+2: ] = param.ωcj # Loss

    σx = (1/ (β * ω**2.0) ) ** 0.5
    σp = (1/β)**0.5 
    #-------- Nuclear DOF ----------
    x = np.zeros(ndof)
    x[0]  = 0.0           # Nuclear DOF
    x[1] = gran(0, σx[1]) # Cavity DOF 
    x[2:param.nbath + 2]  = gran(0, σx[2:param.nbath + 2], param.nbath) # Solvent DOF 

    cQc = x[1] * ccj/ ωcj**2
    x[param.nbath + 2 :]  = gran(- cQc, σx[param.nbath + 2 :], param.nloss) # Loss
    #-------------------------------
    p =  gran(0,σp,ndof)
    return x, p



# Bath -----------------------------------------------
# (1/π) * ∫J(ω)/ω
def FD(ω):
    𝛾  = 0.0009113 # 200 cm-1
    ωb = 0.004556335 
    ɑ  =  0.1 * ωb
    λ  = ɑ * 𝛾 /2 
    return (2 * λ / π) * np.arctan(ω/𝛾)

def discretize2(N, F, ωmax = 0.5/27.2114):
    ω = np.linspace(0.00000001, ωmax, N * 1000)
    dω = ω[1] - ω[0]
    Fω = F(ω)
    λs = Fω[-1] 

    ωj = np.zeros((N))
    cj = np.zeros((N))
    for i in range(N):
        j = i+1
        ωj[i] = ω[np.argmin(np.abs(Fω - ((j-0.5)/N) *  λs))] 
        cj[i] =  ωj[i] * (2 * λs/ N)**0.5 
    return ωj, cj
# Loss -----------------------------------------------
def λc(tau, ωc, g, beta):
    fs =  41.341374575751
    lr = 1.0/(tau*fs)
    return lr/4.0*(1.0-np.exp(-beta*ωc)) * (g*g+ωc*ωc)/(ωc*g)

def discretizeC(N, F, ωc, ωmax = 0.5/27.2114):
    ω = np.linspace(0.00000001, ωmax, N * 1000)
    dω = ω[1] - ω[0]
    Fω = F(ω, ωc)
    λs = Fω[-1] 

    ωj = np.zeros((N))
    cj = np.zeros((N))
    for i in range(N):
        j = i+1
        ωj[i] = ω[np.argmin(np.abs(Fω - ((j-0.5)/N) *  λs))] 
        cj[i] =  ωj[i] * (2 * λs/ N)**0.5 
    return ωj, cj

# (1/π) * ∫J(ω)/ω
def FDc(ω, ωc):
    beta = 1052.8
    𝛾  = 0.0009113  * 5 # 1000 cm-1
    λ  = λc(1000.0, ωc, 𝛾, beta)
    return (2 * λ / π) * np.arctan(ω/𝛾)

# ωcj , ccj = discretizeC(200, FDc, 0.1/27.2114, 0.8/27.2114)  


class param:
    def __init__(self, ωc):
        self.ωc = ωc
        self.nbath = 200
        self.ωbj, self.cbj =  discretize2(self.nbath, FD)# Molecular Bath
        self.nloss = 200
        self.ωcj, self.ccj =  discretizeC(self.nloss, FDc, self.ωc)  

        self.ndof = 1 + self.nbath + self.nloss + 1  

        self.T = 298.0
        self.β = 315774/self.T  
        self.m = np.ones(self.ndof) 

        self.t =  4200 * 3 
        self.dt = 1.0
        self.χ  = self.ωc * 0.0

        self.traj = 1

