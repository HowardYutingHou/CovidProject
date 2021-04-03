import numpy as np
import matplotlib.pyplot as plt 

class regions:

    sir_over_time = []

    onebig_matrix = []

    big_matrix = []
    column_matrix = []
    regions = []

    def __init__(self, *args):
        self.big_matrix = np.array(len(args))
        i = 0
        for region in args:
            # Give each matrix an index, needed for the big_matrix
            region.index = i
            # Add each region to the list of regions
            self.regions.append(region)
            i += 1
        self.big_matrix = [0] * len(self.regions)
        # initialize the two dimensional matrix with numpy
        self.onebig_matrix = np.zeros((3*len(self.regions), 3*len(self.regions)))
        #self.onebig_matrix = [[0] * (3*len(self.regions))] * (3*len(self.regions))
        self.construct_matrix()
        self.construct_column()
        # SIR population at time T0 is added to sir_over_time
        self.sir_over_time.append(self.column_matrix)
        # days to iterate through 
        #self.bruteforce_solver(10)
        self.Henn_solver(10)
        # visualize the solutions from day 0 to day T
        self.plot_solution()
        
    def construct_matrix(self):
        i = 0
        for region in self.regions:
            # initialized a row of the big_matrix
            row_matrix = np.zeros(len(self.regions)).tolist()
            # the part that is on the diagonal of the big_matrix
            diag_matrix = region.SIREpi()
            # store it in the correct position of the row
            row_matrix[region.index] = diag_matrix
            # add the interaction matrices to their proper positions
            for border in list(region.border_InterCoeffs.keys()):
                row_matrix[border.index] = region.border_InterCoeffs[border]
            # add the row to its proper location in the big_matrix
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

    def construct_column(self):
        for region in self.regions:
            self.column_matrix.append([region.S / region.N])
            self.column_matrix.append([region.I / region.N])
            self.column_matrix.append([region.R / region.N])
        self.column_matrix = np.array(self.column_matrix)
    
    
    # Henn Method: second order explicit  
    def Henn_solver(self, Tf):
        
        for t in range(1, Tf + 1):
            A_ut = np.dot(self.onebig_matrix, self.column_matrix)
            A_second = np.dot(self.onebig_matrix, self.onebig_matrix)
            dudt = A_ut + np.dot(A_second,self.column_matrix)/2
            i = 0
            for region in self.regions:
                region.S = dudt[i][0]
                region.I = dudt[i+1][0]
                region.R = dudt[i+2][0]
                i += 3
            # u(i+1) = du/dt + u(i); update the column_matrix
            self.column_matrix = dudt + self.column_matrix
            # update the big_matrix
            self.construct_matrix()
            self.sir_over_time.append(self.column_matrix)
    

    # explicit euler  
    def bruteforce_solver(self, Tf):

        for t in range(1, Tf + 1):
            dudt = np.dot(self.onebig_matrix, self.column_matrix)
            i = 0
            for region in self.regions:
                region.S = dudt[i][0]
                region.I = dudt[i+1][0]
                region.R = dudt[i+2][0]
                i += 3
            # u(i+1) = du/dt + u(i); update the column_matrix
            self.column_matrix = dudt + self.column_matrix
            # update the big_matrix
            self.construct_matrix()
            self.sir_over_time.append(self.column_matrix)
            
    # visualization
    
    def plot_solution(self):
        t = np.arange(11)
        # replace 3 with length of all regions later 
        for i in range(0,3):
            S = []
            I = []
            R = []
            for j in self.sir_over_time:
                S.append(j[i*3])
                I.append(j[i*3+1])
                R.append(j[i*3+2])
            
            # plot every region
            fig = plt.figure(facecolor='w')
            ax = fig.add_subplot(111, axisbelow = True)
            ax.plot(t, S, 'r', alpha=0.5, lw=2, label = 'Susceptible')
            ax.plot(t, I, 'b', alpha=0.5, lw=2, label = 'Infected')
            ax.plot(t, R, 'g', alpha=0.5, lw=2, label = 'Recovered')
            ax.set_xlabel('Time(days)')
            ax.set_ylabel('Number')
            legend = ax.legend()
            legend.get_frame().set_alpha(0.5)
            plt.show()
