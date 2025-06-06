
double
third_body(const double& A, 
           const double& B,
           const double& E, 
           const double& temperature, 
           const double& log_temperature, 
           const double& mixture_concentration) 
{
        return multiply(mixture_concentration,
                        arrhenius(A, B, E, temperature, log_temperature));
}


double
dthird_body_dtemperature(const double& A, 
                         const double& B, 
                         const double& E, 
                         const double& temperature, 
                         const double& log_temperature,
                         const double& mixture_concentration) 
{
        return multiply(mixture_concentration,
                        darrhenius_dtemperature(A,B,E,temperature,log_temperature));
}


double
dthird_body_dlog_temperature(const double& A, 
                             const double& B, 
                             const double& E, 
                             const double& temperature, 
                             const double& log_temperature,
                             const double& mixture_concentration) 
{
        return multiply(mixture_concentration,
                        darrhenius_dlog_temperature(A,B,E,temperature,log_temperature));
}


double
dthird_body_dmixture_concentration(const double& A, 
                                   const double& B, 
                                   const double& E, 
                                   const double& temperature,
                                   const double& log_temperature,
                                   const double& mixture_concentration) 
{
        return arrhenius(A,B,E,temperature, log_temperature);
}