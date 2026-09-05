

Species 
species_specific_heat_constant_pressure_mass_specific(const TemperatureMonomial& temperature_monomial_sequence)  
{
        return
        Species{
        contract(temperature_monomial_sequence, TemperatureMonomial{double(794.9943626577252), double(0.49830932383778), double(-0.0003043787553975861), double(1.1763451247594122e-07), double(-2.6817280987226805e-11), double(3.548428762305633e-15), double(-2.5230195655978164e-19), double(7.391941640834654e-24)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(97.20646377699087), double(5.2426975458232885), double(-0.0033364855500587027), double(1.2747673919011844e-06), double(-3.0170095976914396e-10), double(4.279846635123641e-14), double(-3.2862169066910223e-18), double(1.0310572062221609e-22)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(983.8440000320377), double(0.13804480560710589), double(0.0001521816466804371), double(-1.347143669041383e-07), double(4.365526854598036e-11), double(-7.040180469138235e-15), double(5.690165744581645e-19), double(-1.8589214819088425e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(499.7810613160261), double(1.4509698170607632), double(-0.0010719236358760196), double(4.51929201166496e-07), double(-1.1291461918778277e-10), double(1.6424097823594165e-14), double(-1.2720906002848855e-18), double(4.01984105571123e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(1791.009639026644), double(-0.0011416578621790767), double(0.0009403928864220311), double(-5.741870707526078e-07), double(1.6067943093155351e-10), double(-2.3996473245208053e-14), double(1.8567010154976688e-18), double(-5.862323218475164e-23)})};

}


Species 
species_specific_heat_constant_pressure_mass_specific(const double& temperature)  
{
    return species_specific_heat_constant_pressure_mass_specific(temperature_monomial(temperature));
}


Species 
dspecies_specific_heat_constant_pressure_mass_specific_dtemperature(const double& temperature)  
{
    return species_specific_heat_constant_pressure_mass_specific(dtemperature_monomial_dtemperature(temperature));
}
    

Species 
species_specific_heat_constant_volume_mass_specific(const TemperatureMonomial& temperature_monomial_sequence)  
{
        return
        Species{
        contract(temperature_monomial_sequence, TemperatureMonomial{double(535.1511656406228), double(0.49830932383778), double(-0.0003043787553975861), double(1.1763451247594122e-07), double(-2.6817280987226805e-11), double(3.548428762305633e-15), double(-2.5230195655978164e-19), double(7.391941640834654e-24)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(-56.50319048155457), double(5.2426975458232885), double(-0.0033364855500587027), double(1.2747673919011844e-06), double(-3.0170095976914396e-10), double(4.279846635123641e-14), double(-3.2862169066910223e-18), double(1.0310572062221609e-22)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(687.0049204835464), double(0.13804480560710589), double(0.0001521816466804371), double(-1.347143669041383e-07), double(4.365526854598036e-11), double(-7.040180469138235e-15), double(5.690165744581645e-19), double(-1.8589214819088425e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(310.8546458520701), double(1.4509698170607632), double(-0.0010719236358760196), double(4.51929201166496e-07), double(-1.1291461918778277e-10), double(1.6424097823594165e-14), double(-1.2720906002848855e-18), double(4.01984105571123e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(1329.4796574472246), double(-0.0011416578621790767), double(0.0009403928864220311), double(-5.741870707526078e-07), double(1.6067943093155351e-10), double(-2.3996473245208053e-14), double(1.8567010154976688e-18), double(-5.862323218475164e-23)})};

}


Species 
species_specific_heat_constant_volume_mass_specific(const double& temperature)  
{
    return species_specific_heat_constant_volume_mass_specific(temperature_monomial(temperature));
}


Species 
dspecies_specific_heat_constant_volume_mass_specific_dtemperature(const double& temperature)  
{
    return species_specific_heat_constant_volume_mass_specific(dtemperature_monomial_dtemperature(temperature));
}
    

Species 
species_enthalpy_mass_specific(const TemperatureEnergyMonomial& temperature_energy_monomial_sequence)  
{
        return
        Species{
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-256706.70514570823), double(794.9943626577252), double(0.24915466191889), double(-0.00010145958513252871), double(2.9408628118985305e-08), double(-5.3634561974453606e-12), double(5.914047937176056e-16), double(-3.6043136651397374e-20), double(9.239927051043317e-25)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(1814077.7851878793), double(97.20646377699087), double(2.6213487729116443), double(-0.0011121618500195675), double(3.186918479752961e-07), double(-6.03401919538288e-11), double(7.133077725206069e-15), double(-4.694595580987174e-19), double(1.2888215077777012e-23)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-4246635.755496567), double(983.8440000320377), double(0.06902240280355294), double(5.07272155601457e-05), double(-3.3678591726034574e-08), double(8.731053709196072e-12), double(-1.1733634115230391e-15), double(8.128808206545207e-20), double(-2.323651852386053e-24)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-9146400.682877207), double(499.7810613160261), double(0.7254849085303816), double(-0.00035730787862533987), double(1.12982300291624e-07), double(-2.2582923837556555e-11), double(2.7373496372656944e-15), double(-1.817272286121265e-19), double(5.0248013196390375e-24)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-13964700.873892631), double(1791.009639026644), double(-0.0005708289310895384), double(0.0003134642954740104), double(-1.4354676768815194e-07), double(3.2135886186310705e-11), double(-3.999412207534675e-15), double(2.6524300221395266e-19), double(-7.327904023093956e-24)})};

}


Species 
species_enthalpy_mass_specific(const double& temperature)  
{
    return species_enthalpy_mass_specific(temperature_energy_monomial(temperature));
}


Species 
dspecies_enthalpy_mass_specific_dtemperature(const double& temperature)  
{
    return species_enthalpy_mass_specific(dtemperature_energy_monomial_dtemperature(temperature));
}
    

Species 
species_internal_energy_mass_specific(const TemperatureEnergyMonomial& temperature_energy_monomial_sequence)  
{
        return
        Species{
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-256706.70514570823), double(535.1511656406228), double(0.24915466191889), double(-0.00010145958513252871), double(2.9408628118985305e-08), double(-5.3634561974453606e-12), double(5.914047937176056e-16), double(-3.6043136651397374e-20), double(9.239927051043317e-25)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(1814077.7851878793), double(-56.50319048155457), double(2.6213487729116443), double(-0.0011121618500195675), double(3.186918479752961e-07), double(-6.03401919538288e-11), double(7.133077725206069e-15), double(-4.694595580987174e-19), double(1.2888215077777012e-23)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-4246635.755496567), double(687.0049204835464), double(0.06902240280355294), double(5.07272155601457e-05), double(-3.3678591726034574e-08), double(8.731053709196072e-12), double(-1.1733634115230391e-15), double(8.128808206545207e-20), double(-2.323651852386053e-24)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-9146400.682877207), double(310.8546458520701), double(0.7254849085303816), double(-0.00035730787862533987), double(1.12982300291624e-07), double(-2.2582923837556555e-11), double(2.7373496372656944e-15), double(-1.817272286121265e-19), double(5.0248013196390375e-24)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-13964700.873892631), double(1329.4796574472246), double(-0.0005708289310895384), double(0.0003134642954740104), double(-1.4354676768815194e-07), double(3.2135886186310705e-11), double(-3.999412207534675e-15), double(2.6524300221395266e-19), double(-7.327904023093956e-24)})};

}


Species 
species_internal_energy_mass_specific(const double& temperature)  
{
    return species_internal_energy_mass_specific(temperature_energy_monomial(temperature));
}


Species 
dspecies_internal_energy_mass_specific_dtemperature(const double& temperature)  
{
    return species_internal_energy_mass_specific(dtemperature_energy_monomial_dtemperature(temperature));
}
    

Species 
species_entropy_mass_specific(const TemperatureEnergyMonomial& temperature_entropy_monomial_sequence)  
{
        return
        Species{
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1745.6981956248583), double(0.49830932383778), double(-0.00015218937769879306), double(3.9211504158647075e-08), double(-6.704320246806701e-12), double(7.096857524611266e-16), double(-4.205032609329694e-20), double(1.055991662976379e-24), double(794.9943626577252)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(3443.4395587801964), double(5.2426975458232885), double(-0.0016682427750293514), double(4.249224639670615e-07), double(-7.542523994228599e-11), double(8.559693270247281e-15), double(-5.47702817781837e-19), double(1.4729388660316584e-23), double(97.20646377699087)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1404.2702314607495), double(0.13804480560710589), double(7.609082334021855e-05), double(-4.49047889680461e-08), double(1.091381713649509e-11), double(-1.408036093827647e-15), double(9.483609574302742e-20), double(-2.655602117012632e-24), double(983.8440000320377)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1621.4940109369472), double(1.4509698170607632), double(-0.0005359618179380098), double(1.5064306705549867e-07), double(-2.8228654796945692e-11), double(3.284819564718833e-15), double(-2.1201510004748092e-19), double(5.7426300795874714e-24), double(499.7810613160261)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(240.5684927720049), double(-0.0011416578621790767), double(0.00047019644321101556), double(-1.9139569025086926e-07), double(4.016985773288838e-11), double(-4.7992946490416104e-15), double(3.0945016924961145e-19), double(-8.37474745496452e-24), double(1791.009639026644)})};

}


Species 
species_entropy_mass_specific(const double& temperature)  
{
    return species_entropy_mass_specific(temperature_entropy_monomial(temperature));
}


Species 
dspecies_entropy_mass_specific_dtemperature(const double& temperature)  
{
    return species_entropy_mass_specific(dtemperature_entropy_monomial_dtemperature(temperature));
}
    

Species 
species_gibbs_energy_mole_specific(const TemperatureGibbsMonomial& temperature_gibbs_monomial_sequence)  
{
        return
        Species{
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-8214101.151252372), double(-30420.621247282324), double(-7.972450872080643), double(0.0016232519025353269), double(-3.136724275170972e-07), double(4.290496785146416e-11), double(-3.784754117875191e-15), double(1.922180477619022e-19), double(0.0), double(-25438.229616321893)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(98127095.55638276), double(-181004.4405749134), double(-141.79399782433666), double(0.030079529395629228), double(-5.746226480226571e-06), double(8.159804157916266e-10), double(-7.716848806236937e-14), double(4.232334402779304e-18), double(0.0), double(-5258.09203862499)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-118948267.51145883), double(-11776.138742318219), double(-1.9333175025275178), double(-0.0007104346539198406), double(3.144457847487427e-07), double(-6.113920359864548e-11), double(6.573181831352069e-15), double(-3.794798631088854e-19), double(0.0), double(-27557.470440897374)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-402523947.65274304), double(-49365.46519986712), double(-31.927865339513566), double(0.007862381215211293), double(-1.6574126845113596e-06), double(2.484629737917565e-10), double(-2.4093604037285204e-14), double(1.3329389339985124e-18), double(0.0), double(-21994.864727456992)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-251574086.24317577), double(27931.197249777324), double(0.010283483193578033), double(-0.0028235296414821486), double(8.619983399673523e-07), double(-1.4473199741159678e-10), double(1.4409882183747444e-14), double(-7.963921141473929e-19), double(0.0), double(-32265.038647064994)})};

}


Species 
species_gibbs_energy_mole_specific(const double& temperature)  
{
    return species_gibbs_energy_mole_specific(temperature_gibbs_monomial(temperature));
}


Species 
dspecies_gibbs_energy_mole_specific_dtemperature(const double& temperature)  
{
    return species_gibbs_energy_mole_specific(temperature_gibbs_monomial(temperature));
}
    

Reactions 
gibbs_reaction(const TemperatureMonomial& log_temperature_monomial_sequence)  
{
        return
        Reactions{
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(-166537.84234308443), double(118520.4767789564), double(-36997.61064016703), double(6519.116036840046), double(-696.1158041547768), double(44.8046189872453), double(-1.6011664378201735), double(0.02437528710353279)}),
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(-37779.6412822259), double(27643.45537441301), double(-8948.71550335367), double(1651.4582261102773), double(-186.5507264866835), double(12.833958063628934), double(-0.4957938201263719), double(0.00826780556438974)}),
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(37779.64128222504), double(-27643.45537441224), double(8948.715503353342), double(-1651.4582261102034), double(186.5507264866735), double(-12.833958063628126), double(0.49579382012633466), double(-0.008267805564389049)})};

}


Reactions 
gibbs_reaction(const double& log_temperature)  
{
    return gibbs_reaction(temperature_monomial(log_temperature));
}


Reactions 
dgibbs_reaction_dlog_temperature(const double& log_temperature)  
{
    return gibbs_reaction(dtemperature_monomial_dtemperature(log_temperature)); //functionality is the same
}
    