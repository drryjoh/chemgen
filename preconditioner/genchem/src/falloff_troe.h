
double
f_cent_troe(const double& alpha, 
            const double& T1, 
            const double& T2,
            const double& T3,
            const double& temperature) 
{
        return
        (double(1) - alpha) * exp_gen(-inv_safe_gen(T3) * temperature)
        + alpha * exp_gen(-inv_gen(T1) * temperature)
        + exp_gen(-inv_gen(temperature) * T2);
}


double
df_cent_troe_dtemperature(const double& alpha, 
                          const double& T1, 
                          const double& T2,
                          const double& T3,
                          const double& temperature) 
{
        return
        (double(1) - alpha) * exp_chain(-inv_safe_gen(T3) * temperature, -inv_safe_gen(T3))
        + alpha * exp_chain(-inv_gen(T1) * temperature, -inv_gen(T1))
        + exp_chain(-inv_gen(temperature) * T2, -dinv_da(temperature) * T2);
}


double
f_troe(const double& alpha, 
       const double& T1, 
       const double& T2,
       const double& T3,
       const double& temperature,
       const double& Pr) 
{
        double f_cent_troe_ = f_cent_troe(alpha, T1, T2, T3, temperature); 
        double log_f_cent_troe_ = log10_gen(f_cent_troe_);
        double log_Pr = log10_gen(Pr);
        double c = double(-0.4) + double(-0.67) * log_f_cent_troe_;
        double n = double(0.75) + double(-1.27) * log_f_cent_troe_;
        double d = double(0.14);
        double a = log_Pr + c;
        double b = n - d * (log_Pr + c);
        double Pr_factor = inv_gen(double(1) + pow_gen2(divide(a, b)));
        
        return
        pow_gen(double(10), (multiply(log_f_cent_troe_, Pr_factor)));
}


double
df_troe_dtemperature(const double& alpha, 
                     const double& T1, 
                     const double& T2,
                     const double& T3,
                     const double& temperature,
                     const double& Pr) 
{

        double f_cent_troe_ = f_cent_troe(alpha, T1, T2, T3, temperature); 
        double df_cent_troe_dtemperature_ =  df_cent_troe_dtemperature(alpha, T1, T2, T3, temperature); 

        double log_f_cent_troe_ = log10_gen(f_cent_troe_);
        double dlog_f_cent_troe_dtemperature_ = log10_chain(f_cent_troe_, df_cent_troe_dtemperature_);

        double log_Pr = log10_gen(Pr);
        double dlog_Pr_dtemperature = double(0.0);

        double c = double(-0.4) + double(-0.67) * log_f_cent_troe_;
        double dc_dtemperature = double(-0.67) * dlog_f_cent_troe_dtemperature_;

        double n = double(0.75) + double(-1.27) * log_f_cent_troe_;
        double dn_dtemperature =  double(-1.27) * dlog_f_cent_troe_dtemperature_;
        
        double d = double(0.14);
        double dd_dtemperature = double(0.0);

        double a = log_Pr + c;
        double da_dtemperature = dc_dtemperature;

        double b = n - multiply(d, (log_Pr + c));
        double db_dtemperature = dn_dtemperature - multiply(d, dc_dtemperature);

        double Pr_factor = inv_gen(double(1) + pow_gen2(divide(a, b)));
        double dPr_factor_dtemperature = 
        inv_chain(double(1) + pow_gen2(divide(a, b)), 
                  pow_gen2_chain(divide(a,
                                        b), 
                                 divide_chain(a, 
                                              da_dtemperature,
                                              b,
                                              db_dtemperature)));

        return
        dpow_db_chain(double(10), 
                      multiply(log_f_cent_troe_, 
                               Pr_factor), 
                      multiply_chain(log_f_cent_troe_, 
                                     dlog_f_cent_troe_dtemperature_,
                                     Pr_factor,
                                     dPr_factor_dtemperature));
}


double
df_troe_dPr(const double& alpha, 
            const double& T1, 
            const double& T2,
            const double& T3,
            const double& temperature,
            const double& Pr) 
{
        double f_cent_troe_ = f_cent_troe(alpha, T1, T2, T3, temperature); 
        double log_f_cent_troe_ = log10_gen(f_cent_troe_);

        double log_Pr = log10_gen(Pr);
        double dlog_Pr_dPr = dlog10_da(Pr);

        double c = double(-0.4) + double(-0.67) * log_f_cent_troe_;
        double n = double(0.75) + double(-1.27) * log_f_cent_troe_;
        double d = double(0.14);

        double a = log_Pr + c;
        double da_dPr = dlog_Pr_dPr;

        double b = n - d * (log_Pr + c);
        double db_dPr =  -d * (dlog_Pr_dPr);

        double Pr_factor = inv_gen(double(1) + pow_gen2(divide(a, b)));
        double dPr_factor_dPr = inv_chain(double(1) + pow_gen2(divide(a, b)),
                                            pow_gen2_chain(divide(a, b),
                                                           divide_chain(a,
                                                                        da_dPr,
                                                                        b,
                                                                        db_dPr)));
        
        return
        dpow_db_chain(double(10), 
                      multiply(log_f_cent_troe_, 
                                Pr_factor),
                      multiply_chain(log_f_cent_troe_,
                                     double(0.0), 
                                     Pr_factor,
                                     dPr_factor_dPr));
}


double
falloff_troe(const double& A_low,  //constant
             const double& B_low, //constant
             const double& E_low,  //constant
             const double& A_high, //constant
             const double& B_high,  //constant
             const double& E_high, //constant
             const double& alpha,//constant
             const double& T1,//constant
             const double& T2,//constant
             const double& T3,//constant
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
                 f_troe(alpha, T1, T2, T3, temperature, Pr));
}


double
dfalloff_troe_dtemperature(const double& A_low,  //constant
                          const double& B_low, //constant
                          const double& E_low,  //constant
                          const double& A_high, //constant
                          const double& B_high,  //constant
                          const double& E_high, //constant
                          const double& alpha,//constant
                          const double& T1,//constant
                          const double& T2,//constant
                          const double& T3,//constant
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
                       f_troe(alpha, T1, T2, T3, temperature, Pr), 
                       df_troe_dtemperature(alpha, T1, T2, T3, temperature, Pr) + df_troe_dPr(alpha, T1, T2, T3, temperature, Pr) * dPr_dtemperature);
}


double
dfalloff_troe_dlog_temperature(const double& A_low,  //constant
                               const double& B_low, //constant
                               const double& E_low,  //constant
                               const double& A_high, //constant
                               const double& B_high,  //constant
                               const double& E_high, //constant
                               const double& alpha,//constant
                               const double& T1,//constant
                               const double& T2,//constant
                               const double& T3,//constant
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
                       multiply_chain(k_high, dk_high_dlog_temperature, Pr_inv_one_p_Pr, dPr_inv_one_p_Pr_dlog_temperature),
                       f_troe(alpha, T1, T2, T3, temperature, Pr),
                       df_troe_dPr(alpha, T1, T2, T3, temperature, Pr) * dPr_dlog_temperature);
}


double
dfalloff_troe_dmixture_concentration(const double& A_low,  //constant
                                     const double& B_low, //constant
                                     const double& E_low,  //constant
                                     const double& A_high, //constant
                                     const double& B_high,  //constant
                                     const double& E_high, //constant
                                     const double& alpha,//constant
                                     const double& T1,//constant
                                     const double& T2,//constant
                                     const double& T3,//constant
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
                       f_troe(alpha, T1, T2, T3, temperature, Pr),
                       df_troe_dPr(alpha, T1, T2, T3, temperature, Pr) * dPr_dmixture_concentration);
}