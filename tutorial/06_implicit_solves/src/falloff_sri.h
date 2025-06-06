
double
f_sri(const double& a, 
      const double& b,
      const double& c,
      const double& d,
      const double& e,
      const double& temperature,
      const double& Pr) 
{
    double inner_expr = 
    multiply(a, 
             exp_gen(-divide(b, temperature))) + exp_gen(-divide(temperature, c));
    
    double exponent = 
    inv_gen(1.0 + pow_gen2(log10_gen(Pr)));
    
    return 
    multiply(pow_gen(inner_expr, 
                     exponent),
             pow_gen(temperature, 
                     e));
}


double
df_sri_dtemperature(const double& a, 
                    const double& b,
                    const double& c,
                    const double& d,
                    const double& e,
                    const double& temperature,
                    const double& Pr) 
{
    double inner_expr = 
    multiply(a, 
             exp_gen(-divide(b, temperature)))
    + exp_gen(-divide(temperature, c));
    
    double dinner_expr_dtemperature = 
    multiply(a, 
             exp_chain(-divide(b, temperature),
                       -ddivide_db(b, temperature)))
    + exp_chain(-divide(temperature, c), 
                -ddivide_da(temperature, c));

    double exponent = 
    inv_gen(1.0 + pow_gen2(log10_gen(Pr)));
    
    double dexponent_dtemperature = double(0.0);

    
    return 
    multiply_chain(pow_gen(inner_expr, 
                           exponent),
                   pow_gen_chain(inner_expr,
                                 dinner_expr_dtemperature,
                                 exponent,
                                 dexponent_dtemperature),
                   pow_gen(temperature, 
                           e),
                   dpow_da(temperature,
                           e));

}


double
df_sri_dPr(const double& a, 
           const double& b,
           const double& c,
           const double& d,
           const double& e,
           const double& temperature,
           const double& Pr) 
{
    double inner_expr = 
    multiply(a, 
             exp_gen(-divide(b, temperature))) + exp_gen(-divide(temperature, c));
   double dinner_expr_dPr = double(0.0);
    
    double exponent = 
    inv_gen(1.0 + pow_gen2(log10_gen(Pr)));
    
    double dexponent_dPr = 
     inv_chain(1.0 + pow_gen2(log10_gen(Pr)),
                     pow_gen2_chain(log10_gen(Pr),
                                    dlog10_da(Pr)));

    return 
    multiply(dpow_db(inner_expr, exponent) * dexponent_dPr,
             pow_gen(temperature, 
                     e));

}


double
falloff_sri(const double& A_low,  //constant
            const double& B_low, //constant
            const double& E_low,  //constant
            const double& A_high, //constant
            const double& B_high,  //constant
            const double& E_high, //constant
            const double& a, //constant
            const double& b, //constant
            const double& c, //constant
            const double& d, //constant
            const double& e,// constant
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

    return
        multiply(multiply(k_high, Pr_inv_one_p_Pr),
                 f_sri(a, b, c, d, e, temperature, Pr));

}


double
dfalloff_sri_dtemperature(const double& A_low,  //constant
                          const double& B_low, //constant
                          const double& E_low,  //constant
                          const double& A_high, //constant
                          const double& B_high,  //constant
                          const double& E_high, //constant
                          const double& a, //constant
                          const double& b, //constant
                          const double& c, //constant
                          const double& d, //constant
                          const double& e,// constant
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
        multiply_chain(multiply(k_high, Pr_inv_one_p_Pr),
                       multiply_chain(k_high, dk_high_dtemperature, Pr_inv_one_p_Pr, dPr_inv_one_p_Pr_dtemperature),
                       f_sri(a, b, c, d, e, temperature, Pr),
                       df_sri_dtemperature(a, b, c, d, e, temperature, Pr) + df_sri_dPr(a, b, c, d, e, temperature, Pr) * dPr_dtemperature);
}


double
dfalloff_sri_dlog_temperature(const double& A_low,  //constant
                              const double& B_low, //constant
                              const double& E_low,  //constant
                              const double& A_high, //constant
                              const double& B_high,  //constant
                              const double& E_high, //constant
                              const double& a, //constant
                              const double& b, //constant
                              const double& c, //constant
                              const double& d, //constant
                              const double& e,// constant
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
        multiply_chain(multiply(k_high, Pr_inv_one_p_Pr),
                       multiply_chain(k_high, dk_high_dlog_temperature,  Pr_inv_one_p_Pr, dPr_inv_one_p_Pr_dlog_temperature),
                       f_sri(a, b, c, d, e, temperature, Pr),
                       df_sri_dPr(a, b, c, d, e, temperature, Pr) * dPr_dlog_temperature);
}


double
dfalloff_sri_dmixture_concentration(const double& A_low,  //constant
                                    const double& B_low, //constant
                                    const double& E_low,  //constant
                                    const double& A_high, //constant
                                    const double& B_high,  //constant
                                    const double& E_high, //constant
                                    const double& a, //constant
                                    const double& b, //constant
                                    const double& c, //constant
                                    const double& d, //constant
                                    const double& e,// constant
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
        divide(k_low,
               k_high);

        double inv_one_p_Pr = inv_gen(double(1)+Pr);
        double dinv_one_p_Pr_dmixture_concentration = inv_chain(double(1)+Pr, dPr_dmixture_concentration);

        double Pr_inv_one_p_Pr = multiply(Pr, inv_one_p_Pr);
        double dPr_inv_one_p_Pr_dmixture_concentration = multiply_chain(Pr, dPr_dmixture_concentration, inv_one_p_Pr, dinv_one_p_Pr_dmixture_concentration);

    return
        multiply_chain(multiply(k_high, Pr_inv_one_p_Pr),
                       multiply(k_high, dPr_inv_one_p_Pr_dmixture_concentration),
                       f_sri(a, b, c, d, e, temperature, Pr),
                       df_sri_dPr(a, b, c, d, e, temperature, Pr) * dPr_dmixture_concentration);

}