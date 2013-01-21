double
sim_alpha_func (alpha)
     double alpha;
{
  double answer;
  answer =
    ((alpha + 1.) / (alpha + 2.)) * ((pow (sim_numax, (alpha + 2.)) -
				      pow (sim_numin,
					   (alpha + 2.))) / (pow (sim_numax,
								  (alpha +
								   1.)) -
							     pow (sim_numin,
								  (alpha +
								   1.))));
  answer = answer - sim_meanfreq;
//      printf("NSH alpha=%.3f,f1=%10.2e,f2=%10.2e,meanfreq=%10.2e,ans=%.3f\n",alpha,sim_numin,sim_numax,sim_meanfreq,answer);
  return (answer);
}
