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



class SIREndemic(object):

    def __init__(self):
        pass

    def SIREnd(self, beta, gamma, mu, S, I, R, N):
        dsdt = [mu/S - mu/N - beta * (1/N) * (I/N), 0, 0]
        didt = [0, beta * (1/N) * (S/N) - (gamma + mu) * (1/N), 0]
        drdt = [0, gamma * I / N, -mu/N]
        return [dsdt, didt, drdt]

class SEIR(object):
    pass
