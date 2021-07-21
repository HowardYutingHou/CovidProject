import numpy as np
import matplotlib.pyplot as plt


class regions:
    t = 0
    sir_over_time = []

    onebig_matrix = []

    big_matrix = []
    column_matrix = []
    regions = []

    lockdown_matrix = []

    def __init__(self, Tf, h, *args):
        self.Tf = Tf
        self.h = h
        self.big_matrix = np.array(len(args))
        self.sir_over_time = []
        i = 0
        self.regions = []
        for region in args:
            # Give each matrix an index, needed for the big_matrix
            region.index = i
            # Add each region to the list of regions
            self.regions.append(region)
            i += 1

        # self.lockdown_matrix = np.zeros((len(self.regions), len(self.regions)))
        # Populate the lockdown matrix
        '''
        for region in self.regions:
            for border in list(region.border_lockdown.keys()):
                self.lockdown_matrix[region.index][border.index] = region.border_lockdown[border]
        '''
        self.big_matrix = [0] * len(self.regions)
        # initialize the two dimensional matrix with numpy
        self.onebig_matrix = np.zeros((3 * len(self.regions), 3 * len(self.regions)))
        # self.onebig_matrix = [[0] * (3*len(self.regions))] * (3*len(self.regions))
        self.construct_matrix()
        self.construct_column()
        # SIR population at time T0 is added to sir_over_time
        self.sir_over_time.append(self.column_matrix)

    def construct_matrix(self):
        i = 0
        for region in self.regions:
            # initialized a row of the big_matrix
            if(self.t > region.vaccination):
                region.gamma += 0.005;

            row_matrix = np.zeros(len(self.regions)).tolist()
            diag_matrix = region.SIR()
            # store it in the correct position of the row
            row_matrix[region.index] = diag_matrix

            titf = region.lockdown

            if ((self.t > titf[0]) & (self.t <= titf[1])):

                region.beta = region.beta / 2
                row_matrix[region.index] = region.SIR()
                region.beta = region.beta * 2  # recover it back to what it was
                for border in list(region.border_InterCoeffs.keys()):
                    if border.index != region.index:
                        row_matrix[border.index] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

            else:
                for border in list(region.border_InterCoeffs.keys()):
                    if border.index != region.index:
                        titf = border.lockdown
                        if ((self.t > titf[0]) & (self.t <= titf[1])):
                            row_matrix[border.index] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
                        else:
                            row_matrix[region.index][0][0] -= region.border_out[border][0]
                            row_matrix[region.index][1][1] -= region.border_out[border][1]
                            row_matrix[region.index][2][2] -= region.border_out[border][2]
                            row_matrix[border.index] = region.border_InterCoeffs[border]

            self.big_matrix[region.index] = row_matrix

        # Copying the nested block matrix to a 2D matrix
        i_start = 0
        for row_region in self.big_matrix:
            j_start = 0
            for column_region in row_region:
                i = i_start
                for row in column_region:
                    j = j_start
                    for item in row:
                        self.onebig_matrix[i][j] = item
                        j += 1
                    i += 1
                j_start += 3
            i_start += 3

    '''
    def construct_column(self):
        self.column_matrix = []
        for region in self.regions:
            self.column_matrix.append([region.S / region.N])
            self.column_matrix.append([region.I / region.N])
            self.column_matrix.append([region.R / region.N])
        self.column_matrix = np.array(self.column_matrix)
    '''

    # without normalization
    def construct_column(self):
        self.column_matrix = []
        for region in self.regions:
            self.column_matrix.append([region.S])
            self.column_matrix.append([region.I])
            self.column_matrix.append([region.R])
        self.column_matrix = np.array(self.column_matrix)

    '''
    # Henn Method: second order explicit
    def Heun_solver(self):
        N = int(self.Tf / self.h) + 1
        for t in range(1, N):
            self.t += 1
            # get A(ui)ui
            dudt = np.dot(self.onebig_matrix, self.column_matrix)
            # temporary u_i+1
            un = self.column_matrix + self.h * dudt
            # update the onebig_matrix
            i = 0
            for region in self.regions:
                region.S = region.N * un[i][0]
                region.I = region.N * un[i + 1][0]
                region.R = region.N * un[i + 2][0]
                i += 3
            self.construct_matrix()
            # get A(u_i+1)u_i+1
            dudt_new = np.dot(self.onebig_matrix, un)
            # finally, update the column_matrix
            self.column_matrix = self.column_matrix + ((self.h) / 2) * (dudt + dudt_new)
            self.sir_over_time.append(self.column_matrix)
    # explicit euler
    def bruteforce_solver(self):
        N = int(self.Tf / self.h) + 1
        for t in range(1, N):
            dudt = (self.h) * np.dot(self.onebig_matrix, self.column_matrix)
            # print(dudt)
            i = 0
            # u(i+1) = du/dt + u(i); update the column_matrix
            self.column_matrix = dudt + self.column_matrix
            # update the big_matrix
            for region in self.regions:
                region.S = region.N * self.column_matrix[i][0]
                region.I = region.N * self.column_matrix[i + 1][0]
                region.R = region.N * self.column_matrix[i + 2][0]
                i += 3
            self.construct_matrix()
            self.sir_over_time.append(self.column_matrix)
    '''

    # Henn Method: without normalization
    def Heun_solver(self):
        N = int(self.Tf / self.h) + 1
        for t in range(1, N):
            self.t += 1
            # get A(ui)ui
            dudt = np.dot(self.onebig_matrix, self.column_matrix)
            # temporary u_i+1
            un = self.column_matrix + self.h * dudt
            # update the onebig_matrix
            i = 0
            for region in self.regions:
                region.S = un[i][0]
                region.I = un[i + 1][0]
                region.R = un[i + 2][0]
                i += 3
            self.construct_matrix()
            # get A(u_i+1)u_i+1
            dudt_new = np.dot(self.onebig_matrix, un)
            # finally, update the column_matrix
            self.column_matrix = self.column_matrix + ((self.h) / 2) * (dudt + dudt_new)
            self.sir_over_time.append(self.column_matrix)

    # explicit euler without normalization
    def bruteforce_solver(self):
        N = int(self.Tf / self.h) + 1
        for t in range(1, N):
            dudt = (self.h) * np.dot(self.onebig_matrix, self.column_matrix)
            # print(dudt)
            i = 0
            # u(i+1) = du/dt + u(i); update the column_matrix
            self.column_matrix = dudt + self.column_matrix
            # update the big_matrix
            for region in self.regions:
                region.S = self.column_matrix[i][0]
                region.I = self.column_matrix[i + 1][0]
                region.R = self.column_matrix[i + 2][0]
                i += 3
            self.construct_matrix()
            self.sir_over_time.append(self.column_matrix)

    # visualization

    def plot_solution(self):
        t = np.arange(0, self.Tf + 1, self.h)
        for i in range(0, len(self.regions)):
            S = []
            I = []
            R = []
            for j in self.sir_over_time:
                S.append(j[i * 3])
                I.append(j[i * 3 + 1])
                R.append(j[i * 3 + 2])

            # plot every region
            fig = plt.figure(facecolor='w')
            ax = fig.add_subplot(111, axisbelow=True)
            ax.plot(t, S, 'r', alpha=0.5, lw=2, label='Susceptible')
            ax.plot(t, I, 'b', alpha=0.5, lw=2, label='Infected')
            ax.plot(t, R, 'g', alpha=0.5, lw=2, label='Recovered')
            ax.set_xlabel('Time(days)')
            ax.set_ylabel('Number')
            legend = ax.legend()
            legend.get_frame().set_alpha(0.5)
            name = self.regions[i].name
            plt.title(name)
            '''
            fig1 = plt.gcf()
            plt.draw()
            plt.show()
            fig1.savefig(name+'new.png', dpi=100)
            '''
