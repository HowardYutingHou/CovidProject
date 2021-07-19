from Region import RegionEpi
from Regions import regions
#import numpy as np


def main():
    # Initialize three regions
    Georgia = RegionEpi(30000, 100, 0, 30100, 0.3, 0.1, "Georgia",[20,50])
    Florida = RegionEpi(45000, 300, 0, 45300, 0.3, 0.1, "Florida",[0,0])
    Alabama = RegionEpi(78000, 600, 0, 78600, 0.3, 0.1, "Alabama",[0,0])
    # Set borders
    Georgia.set_borders(Florida, Alabama)
    Florida.set_borders(Georgia, Alabama)
    Alabama.set_borders(Georgia, Florida)
    # Set interaction coefficients between the borders
    #Georgia.set_InterCoeffs(Florida, 0, 0, 0, 0)
    #Alabama.set_InterCoeffs(Georgia, 0, 0, 0, 0)
    Georgia.set_InterCoeffs(Florida, 0.02, 0.03, 0.015, 0.02, 0.02, 0.03)
    Alabama.set_InterCoeffs(Georgia, 0.015, 0.01, 0.03, 0.01, 0.015, 0.01)
    Florida.set_InterCoeffs(Alabama, 0.01, 0.002, 0.005, 0.01, 0.01, 0.002)

    # Initiate the area, composed of Georgia, Florida, and Alabama.
    area = regions(100, 1, Georgia, Florida, Alabama)
    
    print("\nThe big matrix at last iteration:")
    print(area.onebig_matrix)

    # add lockdown
    #Georgia.lockdown(Florida,20,40)

    # days to iterate through
    # self.bruteforce_solver()
    area.Heun_solver()
    # visualize the solutions from day 0 to day T
    area.plot_solution()
    
    
    print("\nColumn matrix of the du/dt, composed of s, i, r of each region at time Tf:")
    print(area.column_matrix)

    # print("\nThe nested block matrix of the du/dt equation:")
    # print("where each row and column combination represents a 3*3 region-region matrix -- sir matrix for the same region, and interaction matrix for different regions")
    # print("In each 3*3 region matrix, each row and column represents s, i, and r components of the sir model.")
    print("\nThe big matrix at last iteration:")
    print(area.onebig_matrix)

    print("\nThe sir_over_time matrix is meant to record the s, i, r values for the regions at each time t.")
    # print(area.sir_over_time)
    print(len(area.sir_over_time) - 1)


##################
if __name__ == '__main__':
    main()
