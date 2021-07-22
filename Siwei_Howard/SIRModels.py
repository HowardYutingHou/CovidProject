class SIREpidemic(object):

    def __init__(self):
        pass

    # mu: death/birth rate
    # gamma: recovery rate
    # beta: contact rate
    def SIREpi(self, beta, gamma, S, I, R, N):
        dsdt = [-beta * I / N, 0, 0]
        didt = [0, beta * S / N - gamma, 0]
        drdt = [0,  gamma, 0]
        return [dsdt, didt, drdt]
    
    def Vaccine(self, beta, gamma, rho, S, I, R, N):
        dsdt = [-beta * (I/N) - rho, 0, 0]
        didt = [0, beta*(S/N) - gamma, 0]
        drdt = [rho, gamma, 0]
        return [dsdt, didt, drdt]
    

class SIREndemic(object):

    def __init__(self):
        pass


    # how would you even express mu*N, also why divides everything by N 
    def SIREnd(self, beta, gamma, mu, S, I, R, N):
        dsdt = [mu*N/S - mu - beta * (I/N), 0, 0]
        didt = [0, beta*(S/N) - (gamma + mu), 0]
        drdt = [0, gamma, -mu]
        return [dsdt, didt, drdt]
    
    # is it better to have a new class or simply a new function 
    # since vaccine might start on a day a lot later than day first 
    def Vaccine(self, beta, gamma, mu, rho, S, I, R, N):
        dsdt = [mu*N/S - mu - beta * (I/N) - rho, 0, 0]
        didt = [0, beta*(S/N) - (gamma + mu), 0]
        drdt = [rho, gamma, -mu]
        return [dsdt, didt, drdt]

class SEIR(object):
    pass

