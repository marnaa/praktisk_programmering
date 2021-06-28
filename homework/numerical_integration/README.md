The values of the integrations are found in out.stdout.txt 
with clebshaw curtis variable change (var change) and without. 
In piComparePlot.png we measure how CC (clebshaw curtis), No CC and GSL approximate 
the integral int(4*sqrt(1-x*x),0..1). We see CC is better, than no CC. In the bottom 
of stdout.txt the errors belonging to these evaluations are listed. 
The other plots are comparisons between gsl and the homemade CC integrator on the integral 
of exp(-x^2) from different infinite intervals it should yield sqrt(pi) or sqrt(pi)/2.
The measured parameters as a function of absolute accuracy is error, value and evaluations.
eval_c.txt gives the number of evaluations of the CC_integrator. 
It is 832, fun and weird that it doesn't change. 
