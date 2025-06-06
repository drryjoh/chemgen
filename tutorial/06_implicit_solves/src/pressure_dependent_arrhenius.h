
double
log_arrhenius(const double& log_A,
              const double& B,
              const double& E,
              const double& log_temperature,
              const double& inv_universal_gas_constant_temperature) 
{
        return log_A + B * log_temperature - E * inv_universal_gas_constant_temperature;
}


double
dlog_arrhenius_dlog_temperature(const double& log_A,
                                const double& B,
                                const double& E,
                                const double& log_temperature,
                                const double& inv_universal_gas_constant_temperature) 
{
        return B;
}


double
dlog_arrhenius_dinv_universal_gas_constant_temperature(const double& log_A,
                                                       const double& B,
                                                       const double& E,
                                                       const double& log_temperature,
                                                       const double& inv_universal_gas_constant_temperature) 
{
        return - E;
}


double
dlog_arrhenius_diugct(const double& log_A,
                      const double& B,
                      const double& E,
                      const double& log_temperature,
                      const double& inv_universal_gas_constant_temperature) 
{
        return - E;
}


double pressure_dependent_arrhenius(const double& log_k1, //f(temperature)
                                               const double& log_k2, //f(temperature)
                                               const double& log_P, //f(pressure)
                                               const double& log_P1,
                                               const double& log_P2) 
{
    return  exp_gen(log_k1 + (log_k2 - log_k1) * (log_P - log_P1) / (log_P2 - log_P1));
}


double dpressure_dependent_arrhenius_dtemperature(const double& log_k1, 
                                                             const double& dlog_k1_dtemperature, 
                                                             const double& log_k2,
                                                             const double& dlog_k2_dtemperature, 
                                                             const double& log_P, 
                                                             const double& log_P1,
                                                             const double& log_P2) 
{
    return  
    exp_chain(log_k1 + (log_k2 - log_k1) * (log_P - log_P1) / (log_P2 - log_P1),
              dlog_k1_dtemperature + (dlog_k2_dtemperature - dlog_k1_dtemperature) * (log_P - log_P1) / (log_P2 - log_P1));
}


double dpressure_dependent_arrhenius_dpressure(const double& log_k1, 
                                                          const double& log_k2,
                                                          const double& log_P,
                                                          const double& dlog_P_dpressure, 
                                                          const double& log_P1,
                                                          const double& log_P2) 
{
    return  
    exp_chain(log_k1 + (log_k2 - log_k1) * (log_P - log_P1) / (log_P2 - log_P1),
              (log_k2 - log_k1) / (log_P2 - log_P1)) * dlog_P_dpressure;
}