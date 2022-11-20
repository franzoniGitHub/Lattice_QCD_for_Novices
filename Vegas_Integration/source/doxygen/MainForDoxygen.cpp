int main()
{
  #include "../SETTINGS.h" //Include the parameters definition

  double res, err; //result and error outputs
  double xl[N_dim-1]; //lower bound vector
  double xu[N_dim-1]; //upper bound vector

  //initialization of the bound vectors
  for(int i=0; i<N_dim-1; i++) {
  xl[i]=-space_bound;
  xu[i]=space_bound;
  }
  size_t dimension=(size_t)(N_dim-1);

  //Definition and initialization of GSL random generation variables
  //These are necessary for the Vegas routine, which uses Montecarlo importance sampling
  const gsl_rng_type *T;
  gsl_rng *my_gsl_rng;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  my_gsl_rng = gsl_rng_alloc(T);

  //Prepare for the output on standard output with a header
  std::cout << std::fixed << std::setprecision(8) << "x\tpath integral\tsigma\texact asymptotic\n";

  //Open an output file and print a header
  std::ofstream out_file("output_file.dat");
  out_file << std::fixed << std::setprecision(8) << "x\tpath integral\tsigma\texact asymptotic\n";

  //Evaluate the path integral for different cyclic boundary conditions contained in my_fp
  double total_integral = 0.0;
  for(double x_now=-x_loop_bound; x_now<=x_loop_bound+1.e-8; x_now+=x_loop_step)
  {
    size_t calls = 500000; //Number of calls of the integration routine
    struct physical_params my_fp={x_now, time_bound, mass}; //Setting the physical parameters

    //Definition of the function with parameters for Montecarlo integration
    gsl_monte_function my_gsl_integrand = { &integrand, dimension,  &my_fp};

	//Definition and initialization of a workspace to maintain the state of integration
    gsl_monte_vegas_state *my_gsl_state = gsl_monte_vegas_alloc(dimension);

	//Vegas routine to integrate over the hypercubic region whose lower and upper limits are defined by xl and xup, respectively
	//Do 10000 iterations and store result and error in res and error
    gsl_monte_vegas_integrate(&my_gsl_integrand, xl, xu, dimension, 10000, my_gsl_rng, my_gsl_state, &res, &err);
    do //Iterate the integration until the required accuracy target is met, namely, when Chisquare/#dofs is sufficiently close to 1
      {
        gsl_monte_vegas_integrate(&my_gsl_integrand, xl, xu, dimension, calls/5, my_gsl_rng, my_gsl_state, &res, &err);
      }
    while (fabs(gsl_monte_vegas_chisq(my_gsl_state) - 1.0) > 0.5);

	//Print all the results for x=x_now
	std::cout << x_now << "\t" << res << "\t" << err << "\t" << asymptotic(x_now, &my_fp) << std::endl;
	out_file << x_now << "\t" << res << "\t" << err << "\t" << asymptotic(x_now, &my_fp) << std::endl;

	//Update a rough estimate of the integral of the path integral over the start-end points x
    total_integral += x_loop_step*res;

	//Free the memory allocated
    gsl_monte_vegas_free(my_gsl_state);
  }
  std::cout << "Estimate of the zero energy level E_0=" << -(log(total_integral))/time_bound << std::endl;

  //Free the memory allocated and close the output file
  gsl_rng_free(my_gsl_rng);
  out_file.close();
  return 0;
}
