class SIREpidemic(object):

    def __init__(self):
        pass

    # mu: death/birth rate
    # gamma: recovery rate
    # beta: contact rate
    def SIREpi(self, beta, gamma, S, I, R, N):
        dsdt = [-beta * I / N, 0, 0]
        didt = [0, beta * S - gamma, 0]
        drdt = [0,  gamma, 0]
        return [dsdt, didt, drdt]



class SIREndemic(object):
    pass



class SEIR(object):
    pass