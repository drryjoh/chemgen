
ChemicalState 
yass(ChemicalState y,  
     const double& dt, 
     double max_norm  = 0.1,
     double min_dt  = 1e-10,
     int max_iter = 20) 
{

    Species y_init = get_species(y);
    Species y_n = y_init;
    double dt_ = dt;
    double time = 0;
    int iteration = 0;
    const SpeciesJacobian I = jacobian_I();

    while(time < dt && iteration < max_iter)
    {

        double temperature_ = temperature(y[0], y_n);
        SpeciesJacobian J = source_jacobian(y_n, temperature_);
        SpeciesJacobian G = I - scale_gen(dt_, J);
        Species rhs1 = scale_gen(dt_, source(y_n, temperature_));

        #ifdef CHEMGEN_DIRECT_SOLVER
        Species dy = invert_jacobian(G) * rhs1;
        #else
        Species dy = gmres_solve(G, rhs1);
        #endif

        double dy_norm  = norm2(dy);

        if (dy_norm > max_norm)
        {
            dt_ *= double(0.5);
            if (dt_ < min_dt)
            {
                    return set_chemical_state(double(0.0), scale_gen(double(0.0), y_n)); //return garbage
            }
        }
        else
        {
            time+=dt_;
            y_n = y_n + dy; 
        }

        iteration++;
    }

    return 
    set_chemical_state(y[0], 
                        y_n); 
}