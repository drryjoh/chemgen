

Species 
species_specific_heat_constant_pressure_mass_specific(const TemperatureMonomial& temperature_monomial_sequence)  
{
        return
        Species{
        contract(temperature_monomial_sequence, TemperatureMonomial{double(20621.18704899117), double(3.304796009061301e-14), double(-6.094239913692169e-17), double(3.7785770483840835e-20), double(-1.0933786759096274e-23), double(1.6185798428874604e-27), double(-1.1882214489292816e-31), double(3.408835558964188e-36)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(14347.174807166406), double(-1.0913992460265658), double(0.0028959961639176924), double(-1.333165986859689e-06), double(3.104312691669342e-10), double(-3.991981056345451e-14), double(2.6988238681716003e-18), double(-7.517813940388394e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(1453.3409760254697), double(-0.3764015218878373), double(0.000366377065170906), double(-1.8236479633037687e-07), double(5.035676444330227e-11), double(-7.687705315839786e-15), double(6.08552957783078e-19), double(-1.956773167037588e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(794.9943626577177), double(0.49830932383782267), double(-0.00030437875539763876), double(1.1763451247597613e-07), double(-2.6817280987239106e-11), double(3.548428762307942e-15), double(-2.523019565599988e-19), double(7.391941640842675e-24)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(1873.0377487424869), double(-0.64695110817973), double(0.0009701393828009771), double(-5.064305446503823e-07), double(1.3686725706942649e-10), double(-2.028450712472819e-14), double(1.5642569220466262e-18), double(-4.925284399870525e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(1791.0096390266199), double(-0.0011416578620771184), double(0.0009403928864218685), double(-5.741870707525361e-07), double(1.6067943093154131e-10), double(-2.3996473245207703e-14), double(1.8567010154977716e-18), double(-5.86232321847591e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(792.6557571451581), double(1.0853746095603412), double(-0.0006449134517275435), double(2.630436519008229e-07), double(-6.639675592413907e-11), double(9.623326364514323e-15), double(-7.333646382161174e-19), double(2.2863074134307396e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(825.6262960656373), double(1.7088049577478517), double(-0.00098938798717666), double(3.6601133674326856e-07), double(-8.589421542335244e-11), double(1.2201767948677145e-14), double(-9.403190170360256e-19), double(2.9590264851566075e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(5193.160985124949), double(1.9204873841895515e-14), double(-3.015932047804993e-17), double(1.8615234784631648e-20), double(-5.633329588441987e-24), double(8.995704061780883e-28), double(-7.316180957521826e-32), double(2.3959949352364117e-36)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(1004.5620761483399), double(0.03890036930091724), double(0.0002495435053136438), double(-1.7956537214537733e-07), double(5.494811529964073e-11), double(-8.641085321239463e-15), double(6.890583381059709e-19), double(-2.2280133527008652e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(1733.7381393358453), double(0.01957549019317404), double(-6.0401257934852445e-05), double(4.8183656170209146e-08), double(-1.423320659411816e-11), double(2.082373333472383e-15), double(-1.5488933721472951e-19), double(4.730691196439385e-24)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(983.8440000320206), double(0.13804480560717036), double(0.00015218164668034234), double(-1.3471436690408724e-07), double(4.36552685459671e-11), double(-7.040180469136468e-15), double(5.690165744580497e-19), double(-1.8589214819085586e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(499.78106131603596), double(1.4509698170607344), double(-0.001071923635875947), double(4.5192920116646544e-07), double(-1.1291461918777963e-10), double(1.6424097823594746e-14), double(-1.2720906002850242e-18), double(4.019841055711971e-23)}),
        contract(temperature_monomial_sequence, TemperatureMonomial{double(1744.1342575697527), double(-0.31281417518926297), double(0.0007459092002522114), double(-4.248560122881658e-07), double(1.188352803905317e-10), double(-1.784966001862209e-14), double(1.3846109033780516e-18), double(-4.3802982541574154e-23)})};

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
species_enthalpy_mass_specific(const TemperatureEnergyMonomial& temperature_energy_monomial_sequence)  
{
        return
        Species{
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(210118843.07296178), double(20621.18704899117), double(1.6523980045306505e-14), double(-2.0314133045640564e-17), double(9.446442620960209e-21), double(-2.1867573518192548e-24), double(2.6976330714791006e-28), double(-1.6974592127561167e-32), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-4252202.37638016), double(14347.174807166406), double(-0.5456996230132829), double(0.0009653320546392308), double(-3.3329149671492224e-07), double(6.208625383338684e-11), double(-6.653301760575752e-15), double(3.8554626688165716e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(15154842.86934987), double(1453.3409760254697), double(-0.18820076094391866), double(0.00012212568839030198), double(-4.559119908259422e-08), double(1.0071352888660455e-11), double(-1.281284219306631e-15), double(8.6936136826154e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-256706.70508824312), double(794.9943626577177), double(0.24915466191891134), double(-0.00010145958513254625), double(2.9408628118994033e-08), double(-5.363456197447821e-12), double(5.914047937179903e-16), double(-3.60431366514284e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(1655887.8115693242), double(1873.0377487424869), double(-0.323475554089865), double(0.00032337979426699234), double(-1.2660763616259557e-07), double(2.7373451413885296e-11), double(-3.3807511874546982e-15), double(2.2346527457808944e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-13964700.87434836), double(1791.0096390266199), double(-0.0005708289310385592), double(0.0003134642954739562), double(-1.4354676768813402e-07), double(3.213588618630826e-11), double(-3.999412207534617e-15), double(2.652430022139674e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(93173.82839498325), double(792.6557571451581), double(0.5426873047801706), double(-0.00021497115057584782), double(6.576091297520572e-08), double(-1.3279351184827815e-11), double(1.6038877274190537e-15), double(-1.0476637688801678e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-4308855.968629003), double(825.6262960656373), double(0.8544024788739258), double(-0.00032979599572555336), double(9.150283418581714e-08), double(-1.717884308467049e-11), double(2.033627991446191e-15), double(-1.3433128814800366e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-1548340.947715004), double(5193.160985124949), double(9.602436920947758e-15), double(-1.0053106826016643e-17), double(4.653808696157912e-21), double(-1.1266659176883975e-24), double(1.4992840102968137e-28), double(-1.0451687082174037e-32), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-303114.78790659737), double(1004.5620761483399), double(0.01945018465045862), double(8.318116843788127e-05), double(-4.489134303634433e-08), double(1.0989623059928146e-11), double(-1.4401808868732438e-15), double(9.843690544371013e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(59150125.682444215), double(1733.7381393358453), double(0.00978774509658702), double(-2.0133752644950814e-05), double(1.2045914042552286e-08), double(-2.846641318823632e-12), double(3.470622222453971e-16), double(-2.212704817353279e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-4246635.755641076), double(983.8440000320206), double(0.06902240280358518), double(5.0727215560114115e-05), double(-3.367859172602181e-08), double(8.73105370919342e-12), double(-1.1733634115227447e-15), double(8.128808206543568e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-9146400.68256471), double(499.78106131603596), double(0.7254849085303672), double(-0.00035730787862531564), double(1.1298230029161636e-07), double(-2.2582923837555925e-11), double(2.737349637265791e-15), double(-1.8172722861214631e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(24801105.35054808), double(1744.1342575697527), double(-0.15640708759463148), double(0.00024863640008407047), double(-1.0621400307204145e-07), double(2.376705607810634e-11), double(-2.974943336437015e-15), double(1.9780155762543595e-19), double(0.0)})};

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
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(210118843.07296178), double(12372.712229394701), double(1.6523980045306505e-14), double(-2.0314133045640564e-17), double(9.446442620960209e-21), double(-2.1867573518192548e-24), double(2.6976330714791006e-28), double(-1.6974592127561167e-32), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-4252202.37638016), double(10222.937397368172), double(-0.5456996230132829), double(0.0009653320546392308), double(-3.3329149671492224e-07), double(6.208625383338684e-11), double(-6.653301760575752e-15), double(3.8554626688165716e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(15154842.86934987), double(933.654581991265), double(-0.18820076094391866), double(0.00012212568839030198), double(-4.559119908259422e-08), double(1.0071352888660455e-11), double(-1.281284219306631e-15), double(8.6936136826154e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-256706.70508824312), double(535.1511656406153), double(0.24915466191891134), double(-0.00010145958513254625), double(2.9408628118994033e-08), double(-5.363456197447821e-12), double(5.914047937179903e-16), double(-3.60431366514284e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(1655887.8115693242), double(1384.1530178579546), double(-0.323475554089865), double(0.00032337979426699234), double(-1.2660763616259557e-07), double(2.7373451413885296e-11), double(-3.3807511874546982e-15), double(2.2346527457808944e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-13964700.87434836), double(1329.4796574472005), double(-0.0005708289310385592), double(0.0003134642954739562), double(-1.4354676768813402e-07), double(3.213588618630826e-11), double(-3.999412207534617e-15), double(2.652430022139674e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(93173.82839498325), double(540.7481458577182), double(0.5426873047801706), double(-0.00021497115057584782), double(6.576091297520572e-08), double(-1.3279351184827815e-11), double(1.6038877274190537e-15), double(-1.0476637688801678e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-4308855.968629003), double(581.1839306233712), double(0.8544024788739258), double(-0.00032979599572555336), double(9.150283418581714e-08), double(-1.717884308467049e-11), double(2.033627991446191e-15), double(-1.3433128814800366e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-1548340.947715004), double(3115.8965910749685), double(9.602436920947758e-15), double(-1.0053106826016643e-17), double(4.653808696157912e-21), double(-1.1266659176883975e-24), double(1.4992840102968137e-28), double(-1.0451687082174037e-32), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-303114.78790659737), double(707.7653809904459), double(0.01945018465045862), double(8.318116843788127e-05), double(-4.489134303634433e-08), double(1.0989623059928146e-11), double(-1.4401808868732438e-15), double(9.843690544371013e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(59150125.682444215), double(1041.5008053792021), double(0.00978774509658702), double(-2.0133752644950814e-05), double(1.2045914042552286e-08), double(-2.846641318823632e-12), double(3.470622222453971e-16), double(-2.212704817353279e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-4246635.755641076), double(687.0049204835293), double(0.06902240280358518), double(5.0727215560114115e-05), double(-3.367859172602181e-08), double(8.73105370919342e-12), double(-1.1733634115227447e-15), double(8.128808206543568e-20), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(-9146400.68256471), double(310.85464585207995), double(0.7254849085303672), double(-0.00035730787862531564), double(1.1298230029161636e-07), double(-2.2582923837555925e-11), double(2.737349637265791e-15), double(-1.8172722861214631e-19), double(0.0)}),
        contract(temperature_energy_monomial_sequence, TemperatureEnergyMonomial{double(24801105.35054808), double(1255.2495266852204), double(-0.15640708759463148), double(0.00024863640008407047), double(-1.0621400307204145e-07), double(2.376705607810634e-11), double(-2.974943336437015e-15), double(1.9780155762543595e-19), double(0.0)})};

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
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(-3684.452240570594), double(3.304796009061301e-14), double(-3.0471199568460845e-17), double(1.2595256827946944e-20), double(-2.7334466897740686e-24), double(3.237159685774921e-28), double(-1.9803690815488027e-32), double(4.869765084234555e-37), double(20621.18704899117)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(-16715.0161252838), double(-1.0913992460265658), double(0.0014479980819588462), double(-4.44388662286563e-07), double(7.760781729173355e-11), double(-7.983962112690902e-15), double(4.498039780286e-19), double(-1.0739734200554848e-23), double(14347.174807166406)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1883.7545081976687), double(-0.3764015218878373), double(0.000183188532585453), double(-6.078826544345896e-08), double(1.2589191110825568e-11), double(-1.537541063167957e-15), double(1.0142549296384633e-19), double(-2.7953902386251257e-24), double(1453.3409760254697)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1745.698195624891), double(0.49830932383782267), double(-0.00015218937769881938), double(3.921150415865871e-08), double(-6.7043202468097764e-12), double(7.096857524615885e-16), double(-4.2050326093333133e-20), double(1.055991662977525e-24), double(794.9943626577177)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(285.7647400407732), double(-0.64695110817973), double(0.00048506969140048854), double(-1.6881018155012742e-07), double(3.421681426735662e-11), double(-4.056901424945638e-15), double(2.6070948700777103e-19), double(-7.036120571243608e-24), double(1873.0377487424869)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(240.5684927721195), double(-0.0011416578620771184), double(0.00047019644321093424), double(-1.9139569025084536e-07), double(4.016985773288533e-11), double(-4.79929464904154e-15), double(3.094501692496286e-19), double(-8.374747454965586e-24), double(1791.0096390266199)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(2127.929405038866), double(1.0853746095603412), double(-0.00032245672586377175), double(8.76812173002743e-08), double(-1.6599188981034768e-11), double(1.9246652729028644e-15), double(-1.2222743970268625e-19), double(3.2661534477581995e-24), double(792.6557571451581)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1722.759825501238), double(1.7088049577478517), double(-0.00049469399358833), double(1.2200377891442286e-07), double(-2.147355385583811e-11), double(2.440353589735429e-15), double(-1.5671983617267093e-19), double(4.2271806930808676e-24), double(825.6262960656373)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1929.2052430908116), double(1.9204873841895515e-14), double(-1.5079660239024965e-17), double(6.2050782615438825e-21), double(-1.4083323971104968e-24), double(1.7991408123561766e-28), double(-1.2193634929203043e-32), double(3.4228499074805883e-37), double(5193.160985124949)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1094.9449700511468), double(0.03890036930091724), double(0.0001247717526568219), double(-5.985512404845911e-08), double(1.3737028824910182e-11), double(-1.7282170642478926e-15), double(1.1484305635099516e-19), double(-3.182876218144093e-24), double(1004.5620761483399)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(3281.299557629827), double(0.01957549019317404), double(-3.0200628967426223e-05), double(1.606121872340305e-08), double(-3.55830164852954e-12), double(4.1647466669447657e-16), double(-2.5814889535788252e-20), double(6.758130280627693e-25), double(1733.7381393358453)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1404.2702314608314), double(0.13804480560717036), double(7.609082334017117e-05), double(-4.490478896802908e-08), double(1.0913817136491775e-11), double(-1.4080360938272936e-15), double(9.483609574300829e-20), double(-2.6556021170122264e-24), double(983.8440000320206)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(1621.4940109368972), double(1.4509698170607344), double(-0.0005359618179379735), double(1.5064306705548848e-07), double(-2.8228654796944907e-11), double(3.284819564718949e-15), double(-2.1201510004750403e-19), double(5.74263007958853e-24), double(499.78106131603596)}),
        contract(temperature_entropy_monomial_sequence, TemperatureEnergyMonomial{double(658.9819983474936), double(-0.31281417518926297), double(0.0003729546001261057), double(-1.4161867076272192e-07), double(2.9708820097632926e-11), double(-3.569932003724418e-15), double(2.3076848389634195e-19), double(-6.2575689345105936e-24), double(1744.1342575697527)})};

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
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(211799793.81754547), double(24500.08440387826), double(-1.6656171885668956e-14), double(1.0238323055002844e-17), double(-3.1740047206426293e-21), double(5.510628526584521e-25), double(-5.43842827210187e-29), double(2.8517314774302755e-33), double(0.0), double(-20786.1565453831)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-8572439.990782404), double(62621.37691981962), double(1.1001304399947782), double(-0.0009730547110763447), double(2.239718857924277e-07), double(-3.129147193202696e-11), double(2.6826112698641445e-15), double(-1.295435456722368e-19), double(0.0), double(-28923.904411247473)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(242462331.0667286), double(-6886.186101223012), double(3.0110239743417546), double(-0.000976944444278221), double(2.431378647074749e-07), double(-4.028289371641965e-11), double(4.0998532449373605e-15), double(-2.318152088469396e-19), double(0.0), double(-23252.00227543149)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-8214101.149413603), double(-30420.621247283612), double(-7.972450872081325), double(0.0016232519025356077), double(-3.136724275171903e-07), double(4.290496785148384e-11), double(-3.784754117877653e-15), double(1.9221804776206764e-19), double(0.0), double(-25438.229616321652)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(28161684.011359498), double(26994.75205899005), double(5.501348748406334), double(-0.00274986008054937), double(7.177386894057542e-07), double(-1.163850720489868e-10), double(1.1499287089008418e-14), double(-6.334123207915946e-19), double(0.0), double(-31854.752992863476)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-251574086.2513857), double(27931.197249774825), double(0.010283483192659645), double(-0.0028235296414816603), double(8.619983399672447e-07), double(-1.4473199741158582e-10), double(1.4409882183747232e-14), double(-7.963921141474371e-19), double(0.0), double(-32265.038647064557)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(3075295.380004817), double(-44072.042022379734), double(-17.91193718157431), double(0.003547668897953217), double(-7.235015645532132e-07), double(1.0957456630160669e-10), double(-1.0587583666238664e-14), double(5.763198392609803e-19), double(0.0), double(-26162.39592033309)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-146561426.91694692), double(-30515.099870222526), double(-29.061645916417717), double(0.005608840499304487), double(-1.0374591339987947e-06), double(1.4608029217049548e-10), double(-1.3834364500210157e-14), double(7.615240725110328e-19), double(0.0), double(-28082.85283437659)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-6197392.574005971), double(13064.315780977322), double(-3.8434733224659344e-14), double(2.011929274401394e-17), double(-6.209114664953017e-21), double(1.1273988138678536e-24), double(-1.2002074356364103e-28), double(6.972323936413994e-33), double(0.0), double(-20786.15654538309)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-8491457.668415418), double(-2531.9863897932332), double(-0.5448774727979477), double(-0.0011651186263094031), double(4.1919536127338333e-07), double(-7.696582510020675e-11), double(8.069045472973416e-15), double(-4.596019115166826e-19), double(0.0), double(-28141.802001219592)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(710452159.5718374), double(-18587.760195129016), double(-0.11756060635510669), double(0.00012091325150925213), double(-4.8227824521698495e-08), double(8.547752220097659e-12), double(-8.337128702778935e-16), double(4.429466260205038e-20), double(0.0), double(-20823.928791562837)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-118948267.51550652), double(-11776.138742320989), double(-1.9333175025284208), double(-0.0007104346539193983), double(3.144457847486235e-07), double(-6.11392035986269e-11), double(6.573181831350418e-15), double(-3.7947986310880883e-19), double(0.0), double(-27557.470440896897)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(-402523947.6389903), double(-49365.465199864484), double(-31.92786533951293), double(0.007862381215210759), double(-1.6574126845112478e-06), double(2.4846297379174964e-10), double(-2.4093604037286053e-14), double(1.3329389339986578e-18), double(0.0), double(-21994.864727457425)}),
        contract(temperature_gibbs_monomial_sequence, TemperatureGibbsMonomial{double(421792398.69677126), double(18455.184472592962), double(2.6600153387218977), double(-0.002114279628114894), double(6.021271834154029e-07), double(-1.0105158068008861e-10), double(1.011897226455687e-14), double(-5.606685150892982e-19), double(0.0), double(-29662.491318488785)})};

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
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(-31830.465930656686), double(23001.71918278712), double(-7315.0909194787755), double(1318.073324939502), double(-144.4300914547916), double(9.57550526882124), double(-0.35417450689285807), double(0.005618509496681922)}),
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(-31830.465930656686), double(23001.71918278712), double(-7315.0909194787755), double(1318.073324939502), double(-144.4300914547916), double(9.57550526882124), double(-0.35417450689285807), double(0.005618509496681922)}),
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(63912.95310323046), double(-47143.30890067597), double(15392.299936430767), double(-2866.3560338819893), double(327.0302413529644), double(-22.761632880450385), double(0.8918196983721065), double(-0.015133871790663108)}),
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(-31830.465930656686), double(23001.71918278712), double(-7315.0909194787755), double(1318.073324939502), double(-144.4300914547916), double(9.57550526882124), double(-0.35417450689285807), double(0.005618509496681922)}),
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(-31830.465930656686), double(23001.71918278712), double(-7315.0909194787755), double(1318.073324939502), double(-144.4300914547916), double(9.57550526882124), double(-0.35417450689285807), double(0.005618509496681922)}),
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(-31830.465930656686), double(23001.71918278712), double(-7315.0909194787755), double(1318.073324939502), double(-144.4300914547916), double(9.57550526882124), double(-0.35417450689285807), double(0.005618509496681922)}),
        contract(log_temperature_monomial_sequence, TemperatureMonomial{double(-31830.465930656686), double(23001.71918278712), double(-7315.0909194787755), double(1318.073324939502), double(-144.4300914547916), double(9.57550526882124), double(-0.35417450689285807), double(0.005618509496681922)})};

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
    