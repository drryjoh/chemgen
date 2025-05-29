    
double call_forward_reaction_0(const double& temperature, const double& log_temperature)  { return arrhenius(double(103990000.00000003), double(0.0), double(64057040.0), temperature, log_temperature);}
    
double call_forward_reaction_1(const double& temperature, const double& log_temperature)  { return arrhenius(double(5039900000000.001), double(1.5), double(64057040000.0), temperature, log_temperature);}
    
double call_forward_reaction_2(const Species& species, const double& temperature, const double& log_temperature, const double& mixture_concentration)  { return third_body(double(5.903e+24), double(-3.32), double(50542720.0), temperature, log_temperature, mixture_concentration + (double(0.5)-double(1))*species[0] + (double(0.75)-double(1))*species[2] + (double(2.0)-double(1))*species[9]);}
    

double
call_forward_reaction_3(const double& temperature, const double& pressure) 
{
        double log_pressure = log_gen(pressure);
        double rate = double(0);
        /**/if (log_pressure < 9.223503358502464) { return arrhenius(double(0.0005550000000000002), double(2.36), double(-757304.0), temperature); }
        else if (9.223503358502464 <= log_pressure && log_pressure < 11.52608845149651)
        {
        double log_k1 = log_gen(arrhenius(double(0.0005550000000000002), double(2.36), double(-757304.0), temperature));
        double log_k2 = log_gen(arrhenius(double(0.17800000000000002), double(1.68), double(8619040.0), temperature)); 
        return pressure_dependent_arrhenius(log_k1, log_k2, log_pressure,  9.223503358502464, 11.52608845149651);
        }
        else if (11.52608845149651 <= log_pressure && log_pressure < 13.828673544490554)
        {
        double log_k1 = log_gen(arrhenius(double(0.17800000000000002), double(1.68), double(8619040.0), temperature));
        double log_k2 = log_gen(arrhenius(double(2370.0000000000005), double(0.56), double(25133288.0), temperature)); 
        return pressure_dependent_arrhenius(log_k1, log_k2, log_pressure,  11.52608845149651, 13.828673544490554);
        }
        else if (13.828673544490554 <= log_pressure && log_pressure < 16.1312586374846)
        {
        double log_k1 = log_gen(arrhenius(double(2370.0000000000005), double(0.56), double(25133288.0), temperature));
        double log_k2 = log_gen(arrhenius(double(27600000.000000004), double(-0.5), double(47931904.0), temperature)); 
        return pressure_dependent_arrhenius(log_k1, log_k2, log_pressure,  13.828673544490554, 16.1312586374846);
        }

        else { return arrhenius(double(27600000.000000004), double(-0.5), double(47931904.0), temperature); }
        return rate;
}
    
    
double
call_forward_reaction_4(const Species& species, const double& temperature, const double& log_temperature, const double& mixture_concentration)  { return falloff_sri(double(2.2900000000000005e+49), double(-11.3), double(401329280.0), double(5440000000000001.0), double(-1.74), double(361330240.0), double(0.138), double(-670.0), double(0.001), double(1.0), double(0.0), temperature, log_temperature, mixture_concentration + (double(2.0)-double(1))*species[1] + (double(6.0)-double(1))*species[5]);}
    
double
call_forward_reaction_5(const Species& species, const double& temperature, const double& log_temperature, const double& mixture_concentration)  { return falloff_lindemann(double(1588800000000.0002), double(-2.1), double(23012000.0), double(12029000.000000002), double(-0.31), double(29049512.0), temperature, log_temperature, mixture_concentration + (double(0.5)-double(1))*species[1] + (double(1.5)-double(1))*species[9]);}
    
double
call_forward_reaction_6(const Species& species, const double& temperature, const double& log_temperature, const double& mixture_concentration)  { return falloff_troe(double(724170000000.0), double(-1.72), double(2196600.0), double(5286300.000000001), double(0.44), double(0.0), double(0.5), double(90000.0), double(0), double(30.0), temperature, log_temperature, mixture_concentration + (double(2.0)-double(1))*species[1] + (double(0.78)-double(1))*species[3] + (double(13.76974842661546)-double(1))*species[5]);}
    //dcall_forward_reaction_0 temperature unused
    //dcall_forward_reaction_1 temperature unused
    //dcall_forward_reaction_2 temperature unused not needed
    
Species dcall_forward_reaction_2_dspecies(const Species& species, const double& temperature, const double& log_temperature, const double& mixture_concentration)  { Species dmixture_concentration_dspecies = {double(0.5), double(1), double(0.75), double(1), double(1), double(1), double(1), double(1), double(1), double(2.0), double(1), double(1), double(1), double(1)};return scale_gen(dthird_body_dmixture_concentration(double(5.903e+24), double(-3.32), double(50542720.0), temperature, log_temperature, mixture_concentration + (double(0.5)-double(1))*species[0] + (double(0.75)-double(1))*species[2] + (double(2.0)-double(1))*species[9]),dmixture_concentration_dspecies);}
    
//dcall_forward_reaction_3_dtemperature is unused
        
//dcall_forward_reaction_3_dtemperature is unused
        

double
dcall_forward_reaction_3_dpressure(const double& temperature, const double& pressure) 
{
        double inv_universal_gas_constant_temperature = inv_gen(universal_gas_constant() * temperature); 
        double log_pressure = log_gen(pressure);
        double dlog_pressure_dpressure = dlog_da(pressure);
        double rate = double(0);
        /**/if (log_pressure < 9.223503358502464) { return double(0); }
        else if (9.223503358502464 <= log_pressure && log_pressure < 11.52608845149651)
        {
        double log_k1 = log_gen(arrhenius(double(0.0005550000000000002), double(2.36), double(-757304.0), temperature));
        double log_k2 = log_gen(arrhenius(double(0.17800000000000002), double(1.68), double(8619040.0), temperature)); 
        return dpressure_dependent_arrhenius_dpressure(log_k1, log_k2, log_pressure, dlog_pressure_dpressure, 9.223503358502464, 11.52608845149651);
        }
        else if (11.52608845149651 <= log_pressure && log_pressure < 13.828673544490554)
        {
        double log_k1 = log_gen(arrhenius(double(0.17800000000000002), double(1.68), double(8619040.0), temperature));
        double log_k2 = log_gen(arrhenius(double(2370.0000000000005), double(0.56), double(25133288.0), temperature)); 
        return dpressure_dependent_arrhenius_dpressure(log_k1, log_k2, log_pressure, dlog_pressure_dpressure, 11.52608845149651, 13.828673544490554);
        }
        else if (13.828673544490554 <= log_pressure && log_pressure < 16.1312586374846)
        {
        double log_k1 = log_gen(arrhenius(double(2370.0000000000005), double(0.56), double(25133288.0), temperature));
        double log_k2 = log_gen(arrhenius(double(27600000.000000004), double(-0.5), double(47931904.0), temperature)); 
        return dpressure_dependent_arrhenius_dpressure(log_k1, log_k2, log_pressure, dlog_pressure_dpressure, 13.828673544490554, 16.1312586374846);
        }

        else { return double(0); }
        return rate;
}

    //dcall_forward_reaction_4 temperature derivative unused
    
Species
dcall_forward_reaction_4_dspecies(const Species& species, const double& temperature, const double& log_temperature, const double& mixture_concentration)  
    { 
        Species dmixture_concentration_dspecies = {double(1), double(2.0), double(1), double(1), double(1), double(6.0), double(1), double(1), double(1), double(1), double(1), double(1), double(1), double(1)};
        return scale_gen(dfalloff_sri_dmixture_concentration(double(2.2900000000000005e+49), double(-11.3), double(401329280.0), double(5440000000000001.0), double(-1.74), double(361330240.0), double(0.138), double(-670.0), double(0.001), double(1.0), double(0.0), temperature, log_temperature, mixture_concentration + (double(2.0)-double(1))*species[1] + (double(6.0)-double(1))*species[5]), dmixture_concentration_dspecies );
    }
    //dcall_forward_reaction_5 temperature derivative unused
    
Species
dcall_forward_reaction_5_dspecies(const Species& species, const double& temperature, const double& log_temperature, const double& mixture_concentration)  { Species dmixture_concentration_dspecies = {double(1), double(0.5), double(1), double(1), double(1), double(1), double(1), double(1), double(1), double(1.5), double(1), double(1), double(1), double(1)};
return scale_gen(dfalloff_lindemann_dmixture_concentration(double(1588800000000.0002), double(-2.1), double(23012000.0), double(12029000.000000002), double(-0.31), double(29049512.0), temperature, log_temperature, mixture_concentration + (double(0.5)-double(1))*species[1] + (double(1.5)-double(1))*species[9]), dmixture_concentration_dspecies);}
    //dcall_forward_reaction_6 temperature derivative unused
    
Species
dcall_forward_reaction_6_dspecies(const Species& species, const double& temperature, const double& log_temperature, const double& mixture_concentration)  { Species dmixture_concentration_dspecies = {double(1), double(2.0), double(1), double(0.78), double(1), double(13.76974842661546), double(1), double(1), double(1), double(1), double(1), double(1), double(1), double(1)};
return scale_gen(dfalloff_troe_dmixture_concentration(double(724170000000.0), double(-1.72), double(2196600.0), double(5286300.000000001), double(0.44), double(0.0), double(0.5), double(90000.0), double(0), double(30.0), temperature, log_temperature, mixture_concentration + (double(2.0)-double(1))*species[1] + (double(0.78)-double(1))*species[3] + (double(13.76974842661546)-double(1))*species[5]), dmixture_concentration_dspecies);}
