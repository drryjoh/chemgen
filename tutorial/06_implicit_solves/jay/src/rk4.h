
double
rk_source_simple(const double&  y_n) 
{
    return double(3) * y_n;
}


double
rk4_simple(const double&  y_n,
           const double& dt) 
{
    double k1  = rk_source_simple(y_n);

    double y_n_2 = y_n +  scale_gen(dt * double(0.5), k1);
    double k2  = rk_source_simple(y_n_2);

    double y_n_3 = y_n +  scale_gen(dt * double(0.5), k2);
    double k3  = rk_source_simple(y_n_3);

    double y_n_4 = y_n +  scale_gen(dt, k3);
    double k4  = rk_source_simple(y_n_4);

    return y_n +
           scale_gen(dt * double(0.16666666666666666), 
                   k1 + 
                   scale_gen(double(2), k2) +
                   scale_gen(double(2), k3) + 
                   k4);
}

/*

template<typename T>
T rk4(const T& y_n, const double& dt) 
{

    T k1  = source(y_n);
    //std::cout << k1 << std::endl;

    T y_n_2 = y_n +  scale_gen(dt * double(0.5), k1);
    T k2  = source(y_n_2);
    //std::cout << k2 << std::endl;

    T y_n_3 = y_n +  scale_gen(dt * double(0.5), k2);
    T k3  = source(y_n_3);
    //std::cout << k3 << std::endl;

    T y_n_4 = y_n +  scale_gen(dt, k3);
    T k4  = source(y_n_4);
    //std::cout << k4 << std::endl;

    return y_n +
           scale_gen(dt * double(0.16666666666666666), 
                   k1 + 
                   scale_gen(double(2), k2) +
                   scale_gen(double(2), k3) + 
                   k4);
}
*/


ChemicalState
rk4(ChemicalState y, const double& dt) 
{
    ChemicalState y_n = y;
    double temperature_ = temperature(y[0], y_n);

    ChemicalState k1  = source(y_n, temperature_);

    ChemicalState y_n_2 = y_n +  scale_gen(dt * double(0.5), k1);
    temperature_ = temperature(y[0], y_n_2);
    ChemicalState k2  = source(y_n_2, temperature_);
    

    ChemicalState y_n_3 = y_n +  scale_gen(dt * double(0.5), k2);
    temperature_ = temperature(y[0], y_n_3);
    ChemicalState k3  = source(y_n_3, temperature_);

    ChemicalState y_n_4 = y_n +  scale_gen(dt, k3);
    temperature_ = temperature(y[0], y_n_4);
    ChemicalState k4  = source(y_n_4, temperature_);

    return 
    y_n +
           scale_gen(dt * double(0.16666666666666666), 
                   k1 + 
                   scale_gen(double(2), k2) +
                   scale_gen(double(2), k3) + 
                   k4);
}