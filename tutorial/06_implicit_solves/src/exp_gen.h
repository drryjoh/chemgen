//exponent
 double exp_gen(const double& a)  {return std::exp(a);}
 double dexp_da(const double& a)  {return exp_gen(a);}
 double exp_chain(const double& a, const double& a_perturbation)  {return a_perturbation * exp_gen(a);}

//ln
 double log_gen(const double& a)  {return std::log(a);}
 double dlog_da(const double& a)  {return inv_gen(a);}
 double log_chain(const double& a, const double& a_perturbation)  {return inv_gen(a) * a_perturbation;}

//log10
 double log10_gen(const double& a)  {return (a > double(0)) ? std::log10(a) : double(0.0);}
 double dlog10_da(const double& a)  {return inv_gen(a * log_gen(10));}
 double log10_chain(const double& a, const double& a_perturbation)  {return dlog10_da(a) * a_perturbation;}