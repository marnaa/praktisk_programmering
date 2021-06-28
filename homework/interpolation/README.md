The answer for part a) is found in the plots intSin(10,100,500).png where the linear interpolator, is used 
to interpolate a period of 2*pi of sin with a different amount of points to demonstrate the increased accuracy.
Also in the sinus plot is shown the integrated values from 0 to x of sin, yielding -cos(x)+1 which is shown.

For part b) we have done a similar plot with the quadratic spline shown in quadratic_plot.png.
Here the differentiated of sin(x) which is cos(x) is also shown. 

Lastly, for part c) and same plot as in part b) is done with the cubic spline, this is in cubic_plot.png.
Here we also compare with gsl, seeing that it produces pretty much the same result.
This project is formatted a bit unusually and has a main for linear, quadratic, and cubic where
the spline functions are defined and the test are made.