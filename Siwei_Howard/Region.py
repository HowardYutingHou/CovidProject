import math

from SIRModels import SIREpidemic, SIREndemic


class RegionEpi(SIREpidemic):
    name = ""
    index = 0

    rho = 0

    # key is of object region, and value is interaction matrix with the region.
    # border_InterCoeffs = {}

    # key is of object region, and value is move out rate between regions
    # border_out = {}

    # Key is a tuple containing the start date and end date of lockdown.
    # border_lockdown = {}

    def __init__(self, S, I, R, N, beta, gamma, rho_0, rho_K, rho_r, name, lockdown, vac_time):
        self.S = S
        self.I = I
        self.R = R
        self.N = N
        self.beta = beta
        self.gamma = gamma

        self.rho_0 = rho_0
        self.rho_K = rho_K
        self.rho_r = rho_r

        self.name = name
        self.lockdown = lockdown
        self.vac_time = vac_time
        self.border_InterCoeffs = {}
        self.border_out = {}
        # self.border_lockdown = {}

        super(RegionEpi, self).__init__()

    def set_borders(self, *borders):
        for border in borders:
            self.border_InterCoeffs[border] = {}  # can initialize to be anything, which is to be changed later.
            # self.border_lockdown[border] = 0  # default 0 if no lockdown

    def __call__(self):
        print(self.name)

    def SIR(self):
        return super(RegionEpi, self).SIREpi(self.beta, self.gamma, self.S, self.I, self.R, self.N)

    def update_rho(self, t):
        numerator = self.rho_0 * math.exp(self.rho_r * (t - self.vac_time))
        denominator = 1 + self.rho_0 * (math.exp(self.rho_r * (t - self.vac_time)) - 1) / self.rho_K
        self.rho = numerator / denominator

    def SIR_vaccine(self, t):
        self.update_rho(t)
        return super(RegionEpi, self).Vaccine(self.beta, self.gamma, self.rho, self.S, self.I, self.R, self.N)

    # s_i = susceptible inflow rate, s_o = susceptible outflow rate
    # i_i = infected inflow rate, i_o = infected outflow rate
    def set_InterCoeffs(self, border, s_i, s_o, i_i, i_o, r_i, r_o):
        out_list = [s_o, i_o, r_o]
        coeff_matrix = [[s_i, 0, 0], [0, i_i, 0], [0, 0, r_i]]
        self.border_InterCoeffs[border] = coeff_matrix
        self.border_out[border] = out_list
        # here, set the other way around
        if (not border.border_InterCoeffs[self]):
            border.set_InterCoeffs(self, s_o, s_i, i_o, i_i, r_o, r_i)
        # edit the diagonal matrix to include the move out rate from the specific region

    # maybe instead of doing city to city, we just lock down entire city?
    # which means no need for parameter border here and no (if not()) whole thing
    '''
    def lockdown(self, border, ti, tf):
        self.border_lockdown[border] = [ti, tf]
        if (not border.border_lockdown[self]):
            border.lockdown(self, ti, tf)
    '''




class RegionEnd(SIREndemic):
    name = ""
    index = 0

    rho = 0

    # key is of object region, and value is interaction matrix with the region.
    # border_InterCoeffs = {}

    # key is of object region, and value is move out rate between regions
    # border_out = {}

    # Key is a tuple containing the start date and end date of lockdown.
    # border_lockdown = {}

    def __init__(self, S, I, R, N, beta, gamma, mu, rho_0, rho_K, rho_r, name, lockdown, vac_time):
        self.S = S
        self.I = I
        self.R = R
        self.N = N
        self.beta = beta
        self.gamma = gamma
        self.mu = mu

        self.rho_0 = rho_0
        self.rho_K = rho_K
        self.rho_r = rho_r

        self.name = name
        self.lockdown = lockdown
        self.vac_time = vac_time
        self.border_InterCoeffs = {}
        self.border_out = {}
        # self.border_lockdown = {}

        super(RegionEnd, self).__init__()

    def set_borders(self, *borders):
        for border in borders:
            self.border_InterCoeffs[border] = {}  # can initialize to be anything, which is to be changed later.
            # self.border_lockdown[border] = 0  # default 0 if no lockdown

    def __call__(self):
        print(self.name)

    def SIR(self):
        return super(RegionEnd, self).SIREnd(self.beta, self.gamma, self.mu, self.S, self.I, self.R, self.N)

    def update_rho(self, t):
        numerator = self.rho_0 * math.exp(self.rho_r * (t - self.vac_time))
        denominator = 1 + self.rho_0 * (math.exp(self.rho_r * (t - self.vac_time)) - 1) / self.rho_K
        self.rho = numerator / denominator

    def SIR_vaccine(self, t):
        self.update_rho(t)
        return super(RegionEnd, self).Vaccine(self.beta, self.gamma, self.mu, self.rho, self.S, self.I, self.R, self.N)

    # s_i = susceptible inflow rate, s_o = susceptible outflow rate
    # i_i = infected inflow rate, i_o = infected outflow rate
    def set_InterCoeffs(self, border, s_i, s_o, i_i, i_o, r_i, r_o):
        out_list = [s_o, i_o, r_o]
        coeff_matrix = [[s_i, 0, 0], [0, i_i, 0], [0, 0, r_i]]
        self.border_InterCoeffs[border] = coeff_matrix
        self.border_out[border] = out_list
        # here, set the other way around
        if (not border.border_InterCoeffs[self]):
            border.set_InterCoeffs(self, s_o, s_i, i_o, i_i, r_o, r_i)
        # edit the diagonal matrix to include the move out rate from the specific region

    '''
    # maybe instead of doing city to city, we just lock down entire city?
    def lockdown(self, border, ti, tf):
        self.border_lockdown[border] = [ti, tf]
        if (not border.border_lockdown[self]):
            border.lockdown(self, ti, tf)
    '''