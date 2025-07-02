
double
falloff_lindemann(const double& A_low,  //constant
                  const double& B_low, //constant
                  const double& E_low,  //constant
                  const double& A_high, //constant
                  const double& B_high,  //constant
                  const double& E_high, //constant
                  const double& temperature, 
                  const double& log_temperature,
                  const double& mixture_concentration) 
{
        double k_low  = arrhenius(A_low, B_low, E_low, temperature, log_temperature);
        double k_high = arrhenius(A_high, B_high, E_high, temperature, log_temperature);
        
        double Pr = 
        divide(multiply(k_low,
                        mixture_concentration),
               k_high);
        double inv_one_p_Pr = inv_gen(double(1)+Pr);
        double Pr_inv_one_p_Pr = multiply(Pr, inv_one_p_Pr);

        double return_value = 
        multiply(k_high,
                 Pr_inv_one_p_Pr);

        return return_value;
}


double
dfalloff_lindemann_dtemperature(const double& A_low,  //constant
                                const double& B_low, //constant
                                const double& E_low,  //constant
                                const double& A_high, //constant
                                const double& B_high,  //constant
                                const double& E_high, //constant
                                const double& temperature, 
                                const double& log_temperature,
                                const double& mixture_concentration) 
{
        double k_low  = arrhenius(A_low, B_low, E_low, temperature, log_temperature);
        double dk_low_dtemperature  = darrhenius_dtemperature(A_low, B_low, E_low, temperature, log_temperature);

        double k_high = arrhenius(A_high, B_high, E_high, temperature, log_temperature);
        double dk_high_dtemperature = darrhenius_dtemperature(A_high, B_high, E_high, temperature, log_temperature);
        
        double Pr = 
        divide(multiply(k_low,
                        mixture_concentration),
               k_high);

        double dPr_dtemperature = 
        divide_chain(multiply(k_low,
                              mixture_concentration),
                     multiply(dk_low_dtemperature,
                              mixture_concentration),
                     k_high,
                     dk_high_dtemperature);

        double inv_one_p_Pr = inv_gen(double(1)+Pr);
        double dinv_one_p_Pr_dtemperature = inv_chain(double(1)+Pr, dPr_dtemperature);

        double Pr_inv_one_p_Pr = multiply(Pr, inv_one_p_Pr);
        double dPr_inv_one_p_Pr_dtemperature = multiply_chain(Pr, dPr_dtemperature, inv_one_p_Pr, dinv_one_p_Pr_dtemperature);

        return
        multiply_chain(k_high, dk_high_dtemperature, Pr_inv_one_p_Pr, dPr_inv_one_p_Pr_dtemperature);
}


double
dfalloff_lindemann_dlog_temperature(const double& A_low,  //constant
                                    const double& B_low, //constant
                                    const double& E_low,  //constant
                                    const double& A_high, //constant
                                    const double& B_high,  //constant
                                    const double& E_high, //constant
                                    const double& temperature, 
                                    const double& log_temperature,
                                    const double& mixture_concentration) 
{
        double k_low  = arrhenius(A_low, B_low, E_low, temperature, log_temperature);
        double dk_low_dlog_temperature  = darrhenius_dlog_temperature(A_low, B_low, E_low, temperature, log_temperature);

        double k_high = arrhenius(A_high, B_high, E_high, temperature, log_temperature);
        double dk_high_dlog_temperature = darrhenius_dlog_temperature(A_high, B_high, E_high, temperature, log_temperature);
        
        double Pr = 
        divide(multiply(k_low,
                        mixture_concentration),
               k_high);
        
        double dPr_dlog_temperature = 
        divide_chain(multiply(k_low,
                              mixture_concentration),
                     multiply(dk_low_dlog_temperature,
                              mixture_concentration),
                     k_high,
                     dk_high_dlog_temperature);

        double inv_one_p_Pr = inv_gen(double(1)+Pr);
        double dinv_one_p_Pr_dlog_temperature = inv_chain(double(1)+Pr, dPr_dlog_temperature);

        double Pr_inv_one_p_Pr = multiply(Pr, inv_one_p_Pr);
        double dPr_inv_one_p_Pr_dlog_temperature = multiply_chain(Pr, dPr_dlog_temperature, inv_one_p_Pr, dinv_one_p_Pr_dlog_temperature);

        return
        multiply_chain(k_high,
                       dk_high_dlog_temperature,
                       Pr_inv_one_p_Pr,
                       dPr_inv_one_p_Pr_dlog_temperature);

}


double
dfalloff_lindemann_dmixture_concentration(const double& A_low,  //constant
                                          const double& B_low, //constant
                                          const double& E_low,  //constant
                                          const double& A_high, //constant
                                          const double& B_high,  //constant
                                          const double& E_high, //constant
                                          const double& temperature, 
                                          const double& log_temperature,
                                          const double& mixture_concentration) 
{
        double k_low  = arrhenius(A_low, B_low, E_low, temperature, log_temperature);
        double k_high = arrhenius(A_high, B_high, E_high, temperature, log_temperature);
        
        double Pr = 
        divide(multiply(k_low,
                        mixture_concentration),
               k_high);
        
        double dPr_dmixture_concentration = 
        divide(k_low,  //d/dM(k_low*M) = k_low
               k_high);

        double inv_one_p_Pr = inv_gen(double(1)+Pr);
        double dinv_one_p_Pr_dmixture_concentration = inv_chain(double(1)+Pr, dPr_dmixture_concentration);

        double Pr_inv_one_p_Pr = multiply(Pr, inv_one_p_Pr);
        double dPr_inv_one_p_P_dmixture_concentration = multiply_chain(Pr, dPr_dmixture_concentration, inv_one_p_Pr, dinv_one_p_Pr_dmixture_concentration);

        return
        multiply(k_high,
                 dPr_inv_one_p_P_dmixture_concentration); //d/dM(k_high * Pr_inv_one_p_Pr) = k_high * d(Pr_inv_one_p_P)/dM
}