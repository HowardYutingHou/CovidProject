from Region import RegionEpi
from Regions import regions
import numpy as np

def main():
    # Initialize three regions
    Georgia = RegionEpi(30000, 100, 0, 31000, 0.3, 0.9, "Georgia")
    Florida = RegionEpi(45000, 300, 0, 45030, 0.3, 0.9, "Florida")
    Alabama = RegionEpi(78000, 600, 0, 78060, 0.3, 0.9, "Alabama")
    # Set borders
    Georgia.set_borders(Florida, Alabama)
    Florida.set_borders(Georgia, Alabama)
    Alabama.set_borders(Georgia, Florida)
    # Set interaction coefficients between the borders
    Georgia.set_InterCoeffs(Florida, 0.02, 0.03, 0.015, 0.02)
    Georgia.set_InterCoeffs(Alabama, 0.015, 0.01, 0.032, 0.01)
    Florida.set_InterCoeffs(Alabama, 0.009, 0.002, 0.005, 0.012)
    # Initiate the area, composed of Georgia, Florida, and Alabama.
    area = regions(Georgia, Florida, Alabama)

    print("\nColumn matrix of the du/dt, composed of s, i, r of each region:")
    print(area.column_matrix)

    print("\nThe nested block matrix of the du/dt equation:")
    print("where each row and column combination represents a 3*3 region-region matrix -- sir matrix for the same region, and interaction matrix for different regions")
    print("In each 3*3 region matrix, each row and column represents s, i, and r components of the sir model.")
    print(area.big_matrix)

    print("\nThe sir_over_time matrix is meant to record the s, i, r values for the regions at each time t.")
    print("We haven't solved the differential equation, and we only has data for t=0")
    print(area.sir_over_time)

    # onebig_matrix is meant to transform the nested block matrix to a 2D matrix, used for multiplication in solving the DE.
    # However, after transforming, the onebig_matrix is made up of 0's.
    # I tried to debug, as shown in the beginning of the printed result:
    # I have all the items from the nested block matrix put in the correct position of onebig_matrix, though. I am confused.


##################
if __name__ == '__main__':
    main()
