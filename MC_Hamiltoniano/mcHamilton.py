import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

data = np.genfromtxt('datos_observacionales.dat', delimiter=' ')

def modelo(v, t, sigma, beta, rho):
    x, y, z = v
    dvdt = [sigma*(y-x), x*(rho-z)-y, x*y-beta*z]
    return dvdt

def model(datos,param):
    v=datos[:,1:]
    t=datos[:,0]
    sigma, beta, rho = param
    sol = odeint(modelo, v[0,:], t, args=(sigma, beta, rho))
    return sol
    
def loglikelihood(datos, param):
    """Logaritmo natural de la verosimilitud construida con los datos observacionales y los 
        parametros que describen el modelo.
    """
    d = datos[:,1:] -  model(datos, param)
    d = np.linalg.norm(d, axis=1)
    d = -0.5 * np.sum(d**2)
    return d

def logprior(param):
    """Logaritmo natural de los prior para los parametros.
        Todos corresponden a gaussianas con sigma=10.0.
    """
    d = -0.5 * np.sum(param**2/(10.0)**2)
    return d

def divergence_loglikelihood(datos, param):
    """Divergencia del logaritmo de la funcion de verosimilitud.
    """
    n_param = len(param)
    div = np.ones(n_param)
    delta = 1E-5
    for i in range(n_param):
        delta_parameter = np.zeros(n_param)
        delta_parameter[i] = delta
        
        div[i] = loglikelihood(datos, param + delta_parameter) 
        div[i] = div[i] - loglikelihood(datos, param - delta_parameter)
        div[i] = div[i]/(2.0 * delta)
    return div

def hamiltonian(datos, param, param_momentum):
    """Hamiltoniano: energia cinetica + potencial: K+V
    """
    m = 1000.0                                                    #MASA, parametro por modificar
    K = 0.5 * np.sum(param_momentum**2)/m
    V = -loglikelihood(datos, param)     
    return K + V


def leapfrog_proposal(datos, param, param_momentum):
    """Integracion tipo leapfrog. 
        `param` representa las posiciones (i.e. los parametros).
        `param_momemtum` representa el momentum asociado a los parametros.
    """
    N_steps = 3
    delta_t = 1E-2
    m = 100.0
    new_param = param.copy()
    new_param_momentum = param_momentum.copy()
    for i in range(N_steps):
        new_param_momentum = new_param_momentum + divergence_loglikelihood(datos, param) * 0.5 * delta_t 
        new_param = new_param + (new_param_momentum/m) * delta_t
        new_param_momentum = new_param_momentum + divergence_loglikelihood(datos, param) * 0.5 * delta_t
    new_param_momentum = -new_param_momentum
    return new_param, new_param_momentum


def monte_carlo(datos, N=5000):
    param = [np.random.random(3)]
    param_momentum = [np.random.normal(size=3)]
    for i in range(1,N):
        propuesta_param, propuesta_param_momentum = leapfrog_proposal(datos, param[i-1], param_momentum[i-1])
        energy_new = hamiltonian(datos, propuesta_param, propuesta_param_momentum)
        energy_old = hamiltonian(datos, param[i-1], param_momentum[i-1])
        #En este caso el criterio se hace con la energía.
        r = min(1,np.exp(-(energy_new - energy_old)))
        alpha = np.random.random()
        if(alpha<r):
            param.append(propuesta_param)
        else:
            param.append(param[i-1])
        param_momentum.append(np.random.normal(size=3))    

    param = np.array(param)
    return param

param_chain = monte_carlo(data, N=1000)

n_param  = len(param_chain[0])
bestm = []
bestd=[]
nom=['sigma','beta','rho']
plt.figure(figsize=(15,5))
for i in range(n_param):
    plt.subplot(1,n_param,i+1)
    meani=np.mean(param_chain[:,i])
    bestm.append(meani)
    desvi=np.std(param_chain[:,i])
    bestd.append(desvi)
    lab1='Media: {:.3f} '.format(meani)
    lab2='Desviacion: {:.3f}'.format(desvi)
    _=plt.hist(param_chain[:,i], bins=60,label=lab1+lab2)
    plt.title(nom[i])
    plt.legend()
plt.savefig('histParam.pdf')

resul=model(data,bestm)

plt.figure(figsize=(15,5))
plt.subplot(1,3,1)
plt.plot(data[:,0], data[:,1])
plt.plot(data[:,0], resul[:,0])
plt.xlabel('t')
plt.ylabel('x')
plt.title('Desarrollo de x')

plt.subplot(1,3,2)
plt.plot(data[:,0], data[:,2])
plt.plot(data[:,0], resul[:,1])
plt.xlabel('t')
plt.ylabel('y')
plt.title('Desarrollo de y')

plt.subplot(1,3,3)
plt.plot(data[:,0], data[:,3])
plt.plot(data[:,0], resul[:,2])
plt.xlabel('t')
plt.ylabel('z')
plt.title('Desarrollo de z')
plt.savefig('desarrolloTemporal.pdf')
