In these exercises the gram-schmidt orthogonalisation is introduced. As test a random square matrix is generated.
which is factorized and the nescessary tests are made from ex A) also the inverse is calculated, and it is 
tested whether it is the inverse. All this done in out.gstest.txt. Here, there is a tall test for tall matrices
and a square test, where the equation Ax=b is also solved, and an inverse is found. 
For part c) we test the time of the gsl_routine and see whether it goes as O(N^3) this plot is in timeplot.png.
The gsl time is also tested and plotted in gslplot.png together with a O(N^3) graph which is clearly grows slower than.
finally, in compareplot.png the gsl time and gs_routine time is plotted together. 