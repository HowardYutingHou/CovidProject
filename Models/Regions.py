import numpy as np

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
        self.onebig_matrix = [[0] * (3*len(self.regions))] * (3*len(self.regions))
        self.construct_matrix()
        self.construct_column()
        # SIR population at time T0 is added to sir_over_time
        self.sir_over_time.append(self.column_matrix)


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

        print("DEBUGING DEBUGING DEBUGING")
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
                        print(self.onebig_matrix[i][j], end=", ")
                        print("at i={} and j={}".format(i, j))
                        j += 1
                    i += 1

                j_start += 3

            i_start += 3


        print("However, the onebig_matrix is all made up of 0's!?")
        print(self.onebig_matrix)

    def construct_column(self):
        for region in self.regions:
            self.column_matrix.append([region.S / region.N])
            self.column_matrix.append([region.I / region.N])
            self.column_matrix.append([region.R / region.N])
        self.column_matrix = np.array(self.column_matrix)



    # Not done yet. As I failed to transform the nested block matrix to a 2D matrix, multiplications can't be performed.
    def bruteforce_solver(self, Tf):
        for t in range(1, Tf + 1):
            dudt = np.dot(self.big_matrix, self.column_matrix)
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
