
ChemicalState 
yass(ChemicalState y,  
     const double& dt, 
     int &n_gmres_iter_total,
     double max_norm  = 0.1,
     double min_dt  = 1e-10,
     int max_iter = 20)
{

    ChemicalState y_init = y;
    ChemicalState y_n = y_init;
    double dt_ = dt;
    double time = 0;
    int iteration = 0;

    while(time < dt && iteration < max_iter)
    {

        double temperature_ = temperature(y[0], y_n);
        SparseMatrix<double> G = source_jacobian(y_n, temperature_, -dt_, 1);
        ChemicalState rhs1 = scale_gen(dt_, source(y_n, temperature_));

        #ifdef CHEMGEN_DIRECT_SOLVER
        ChemicalState dy = invert_jacobian(G) * rhs1;
        #else
        ChemicalState dy = gmres_solve(G, rhs1, n_gmres_iter_total);
        #endif

        double dy_norm  = norm2(get_species(dy));

        if (dy_norm > max_norm)
        {
            std::cout << "refining YASS; dy: " <<dy_norm<<std::endl;
            dt_ *= double(0.5);
            if (dt_ < min_dt)
            {
                    return scale_gen(double(0.0), y_n); //return garbage
            }
        }
        else
        {
            time+=dt_;
            y_n = y_n + dy; 
        }

        iteration++;
    }

    return y_n;
}