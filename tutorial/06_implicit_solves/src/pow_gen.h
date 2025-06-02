//pow
 double pow_gen(const double& a, const double& b)  { return std::pow(a, b);}
// Partial derivative of a^b with respect to a
 double dpow_da(const double& a, const double& b)  { return b * std::pow(a, b - double(1));}
//chain rule if b = constant
 double dpow_da_chain(const double& a, const double& a_chain, const double& b)  { return a_chain * b * std::pow(a, b - double(1));}
// Partial derivative of a^b with respect to b
 double dpow_db(const double& a, const double& b)  { return std::pow(a, b) * std::log(a);}
// chain rule if a = constant
 double dpow_db_chain(const double& a, const double& b, const double& b_perturbation)  { return std::pow(a, b) * std::log(a) * b_perturbation;}

 
double
pow_gen_chain(const double& a, 
              const double& a_perturbation,
              const double& b,
              const double& b_perturbation) 
{ 
    return dpow_da(a, b) * a_perturbation + dpow_db(a, b) * b_perturbation;
}

//simple pow
 double pow_gen2(const double& a)  { return a * a;}
 double dpow_gen2_da(const double& a)  { return double(2) * a;}
 double pow_gen2_chain(const double& a, const double& a_perturbation)  { return double(2) * a * a_perturbation;}

 double pow_gen3(const double& a)  { return a * a * a;}
 double dpow_gen3_da(const double& a)  { return double(3) * a * a;}
 double pow_gen3_chain(const double& a, const double& a_perturbation)  { return double(3) * a * a * a_perturbation;}

 double pow_gen4(const double& a)  { return a * a * a * a;}
 double dpow_gen4_da(const double& a)  { return double(4) * a * a * a;}
 double pow_gen4_chain(const double& a, const double& a_perturbation)  { return double(4) * a * a * a * a_perturbation;}