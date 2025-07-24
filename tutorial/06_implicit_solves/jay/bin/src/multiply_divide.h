//scalar multiply/divides
//multiply
 double multiply(const double& a, const double& b)  {return a * b;}
//derivatives
 double dmultiply_da(const double& a, const double& b)  {return b;}
 double dmultiply_db(const double& a, const double& b)  {return a;}
//derivatives * increment
 double multiply_chain(const double& a, const double& a_perturbation, const double& b, const double& b_perturbation)  {return a_perturbation * b + a * b_perturbation;}

//pow
 double pow2(const double& a)  {return a * a;}
 double dpow2_da(const double& a)  {return double(2) * a;}

//iverse
 double inv_gen(const double& a)  {return double(1)/a;}
 double dinv_da(const double& a)  {return double(-1)/pow2(a);}
 double inv_safe_gen(const double& a)  {return (a != double(0)) ? (double(1.0) / a) : double(0.0);}
 double dinv_safe_da(const double& a)  {return (a != double(0)) ? (double(-1.0) / pow2(a)) : double(0.0);}
 double inv_chain(const double& a, const double& a_perturbation)  {return double(-a_perturbation)/pow2(a);}

//divide
 double divide(const double& a, const double& b)  {return a / b;}
//derivatives
 double ddivide_da(const double& a, const double& b)  {return inv_gen(b);}
 double ddivide_db(const double& a, const double& b)  {return -a / pow2(b);}
//derivatives * increment
 double divide_chain(const double& a, const double& a_perturbation, const double& b, const double& b_perturbation)  {return a_perturbation / b + -a / pow2(b) * b_perturbation;}

 double sqrt_gen(const double& a)  {return std::sqrt(a);}
 double dsqrt_da(const double& a)  {return double(0.5) * inv_gen(std::sqrt(a));}
 double sqrt_chain(const double& a, const double& da)  { return double(0.5) * da * inv_gen(std::sqrt(a));}

//abs/min/max
 double abs_gen(const double& a)  { return std::abs(a);}