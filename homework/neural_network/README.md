In this we use the simplex downhill method since it proved to be extremely! robust so the initial parameters
actually has little importance.
It, however, need to be modified to include a cost function, which is done in minimization.c.
It also draws upon the numerical integration previously developed, which is found in
integrator_homemade.c.
In comparePlot.png we train the network to exp(-x^2), and show the 
routine's ability to do find the integral and differential of the function.

In diffepPlot.png we solve the differential equation y''=-y through the neural network,
which is the requirement of c). This add's the requirement of the second derivative to the
ann structure, if not using untrained, this is not necessary to give the ann. 

For the first time the functions are not all in main, primarily because there are so many.
And dimitri told me so :)




