# 2D-Magnetic-Interface-Inversion
This is a matlab program to invert the upper surface of the magnetic bedrock. 
# magnetic_froward_2D.m
This is a magnetic forward function refered to Liu(2020), a analytical solution derived by Guan(2005) is used to verify the correctness of this program in `test_forward.m`.
# Data_produce.m
A program to calculate the magnetic responce caused by a synthetic model.
# Inversion_main.m
Using Gauss-Newton iterative algorithms, this program aims at recovering the magnetic interface based on the observation data in `Data_produce`.
