
double
arrhenius(const double& A, const double& B, const double& E, const double& temperature) 
{
        double pow_term = pow_gen(temperature, B);
        double exp_term = exp_gen(divide(-E, universal_gas_constant() * temperature));

        return multiply(A,
                        multiply(pow_term,
                                 exp_term));
}


double
darrhenius_dtemperature(const double& A, const double& B, const double& E, const double& temperature) 
{
        double pow_term = pow_gen(temperature, B);
        double dpow_term_dtemperature = dpow_da(temperature, B);

        double exp_term =
            exp_gen(divide(-E,
                           universal_gas_constant() * temperature));
        double dexp_term_dtemperature =
            exp_chain(divide(-E,
                             universal_gas_constant() * temperature),
                      ddivide_db(-E,
                                 universal_gas_constant() * temperature) * universal_gas_constant());

        return multiply(A,
                        multiply_chain(pow_term,
                                       dpow_term_dtemperature,
                                       exp_term,
                                       dexp_term_dtemperature));
}


double
arrhenius(const double& A, const double& B, const double& E, const double& temperature, const double& log_temperature) 
{

        double exp_term = exp_gen(B * log_temperature + divide(-E, universal_gas_constant() * temperature));

        return multiply(A, exp_term);
}


double
darrhenius_dtemperature(const double& A, const double& B, const double& E, const double& temperature, const double& log_temperature) 
{

        double dexp_term_dtemperature = exp_chain(B * log_temperature + divide(-E, universal_gas_constant() * temperature),
                                                    ddivide_db(-E, universal_gas_constant() * temperature) * universal_gas_constant());

        return multiply(A, dexp_term_dtemperature);
}


double
darrhenius_dlog_temperature(const double& A, const double& B, const double& E, const double& temperature, const double& log_temperature) 
{

        double dexp_term_dlog_temperature = exp_chain(B * log_temperature + divide(-E, universal_gas_constant() * temperature),
                                                    B);

        return multiply(A, dexp_term_dlog_temperature);
}