#pragma once
#include <numeric>  // For std::accumulate

double 
sum_gen(const Species& a_s)  {return std::accumulate(a_s.begin(), a_s.end(), double(0.0));}


Species 
fill_gen(const double& a) 
{
    Species b_s{};  // Create an empty array
    std::fill(b_s.begin(), b_s.end(), a);  // Fill it with 'value'
    return b_s;  // Return the filled array
}


Species 
scale_gen(const double& a, const Species& b_s) 
{
    Species c_s;
    std::transform(b_s.begin(), b_s.end(), c_s.begin(), [&a](double s) { return a * s; });
    return c_s;
}


ChemicalState 
scale_gen(const double& a, const ChemicalState& b_s) 
{
    ChemicalState c_s;
    std::transform(b_s.begin(), b_s.end(), c_s.begin(), [&a](double s) { return a * s; });
    return c_s;
}


Reactions 
scale_gen(const double& a, const Reactions& b_s) 
{
    Reactions c_s;
    std::transform(b_s.begin(), b_s.end(), c_s.begin(), [&a](double s) { return a * s; });
    return c_s;
}


TemperatureMonomial 
scale_gen(const double& a, const TemperatureMonomial& b_s) 
{
    TemperatureMonomial c_s;
    std::transform(b_s.begin(), b_s.end(), c_s.begin(), [&a](double s) { return a * s; });
    return c_s;
}


double 
scale_gen(const double& a, const double& b) 
{
    return a * b;
}


SpeciesJacobian 
scale_gen(const double& a, const SpeciesJacobian& b) 
{
    SpeciesJacobian c = {};

    for (int i = 0; i<n_species; i++)
    {
        for (int j = 0; j<n_species; j++)
        {
            c[i][j] = a * b[i][j];
        }
    }
    return c;
}


Species 
operator+(const Species& a_s, const Species& b_s) 
{
    Species c_s;
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a + b; });
    return c_s;
}


SpeciesJacobian 
operator+(const SpeciesJacobian& a, const SpeciesJacobian& b) 
{
    SpeciesJacobian c = {};

    for (int i = 0; i<n_species; i++)
    {
        for (int j = 0; j<n_species; j++)
        {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    return c;
}


Species 
operator-(const Species& a_s, const Species& b_s_positive) 
{
    Species c_s;
    Species b_s = scale_gen(double(-1.0), b_s_positive);
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a + b; });
    return c_s;
}


Species 
safe_inv(const Species& a_s) 
{
    Species b = {};
    for (int i = 0; i<n_species; i++)
    {
        b[i] = inv_safe_gen(a_s[i]);
    }
    return b;
}


Species 
safe_divide(const Species& a_s, const Species& b_s_divide) 
{
    Species c_s;
    Species b_s = safe_inv(b_s_divide);
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a * b; });
    return c_s;
}


Reactions 
operator+(const Reactions& a_s, const Reactions& b_s) 
{
    Reactions c_s;
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a + b; });
    return c_s;
}


Reactions 
operator-(const Reactions& a_s, const Reactions& b_s_positive) 
{
    Reactions c_s;
    Reactions b_s = scale_gen(double(-1.0), b_s_positive);
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a + b; });
    return c_s;
}


SpeciesJacobian 
operator-(const SpeciesJacobian& a, const SpeciesJacobian& b) 
{
    SpeciesJacobian c = {};

    for (int i = 0; i<n_species; i++)
    {
        for (int j = 0; j<n_species; j++)
        {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    return c;
}


double 
L2(const SpeciesJacobian& a) 
{
    double sum = double(0.0);

    for (int i = 0; i<n_species; i++)
    {
        for (int j = 0; j<n_species; j++)
        {
            sum+= pow2(a[i][j]);
        }
    }
    return sqrt_gen(sum);
}



TemperatureMonomial 
operator+(const TemperatureMonomial& a_s, const TemperatureMonomial& b_s) 
{
    TemperatureMonomial c_s;
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a + b; });
    return c_s;
}


TemperatureMonomial 
operator-(const TemperatureMonomial& a_s, const TemperatureMonomial& b_s_positive) 
{
    TemperatureMonomial c_s;
    TemperatureMonomial b_s = scale_gen(double(-1.0), b_s_positive);
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a + b; });
    return c_s;
}



ChemicalState 
operator+(const ChemicalState& a_s, const ChemicalState& b_s) 
{
    ChemicalState c_s;
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a + b; });
    return c_s;
}


ChemicalState 
operator-(const ChemicalState& a_s, const ChemicalState& b_s_positive) 
{
    ChemicalState c_s;
    ChemicalState b_s = scale_gen(double(-1.0), b_s_positive);
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a + b; });
    return c_s;
}


Species 
operator*(const Species& a_s, const Species& b_s) 
{
    Species c_s;
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a * b; });
    return c_s;
}


ChemicalState 
operator*(const ChemicalState& a_s, const ChemicalState& b_s) 
{
    ChemicalState c_s;
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a * b; });
    return c_s;
}


Species 
get_species(const ChemicalState& y) 
{
    Species s;
    for (int i = 1; i <= n_species; i++)
    {
        s[i-1] = y[i];
    }
    return s;
}


double 
get_energy(const ChemicalState& y) 
{
    return y[0];
}


ChemicalState 
set_chemical_state(const double& internal_energy,
                   const Species& species) 
{
    ChemicalState y;
    for (int i = 1; i <= n_species; i++)
    {
        y[i] = species[i-1];
    }
    y[0] = internal_energy;
    return y;
}


double norm2(const Species& x) 
{
    double sum = double(0.0);
    for (int i  = 0; i<n_species; i++)
    {
        sum += x[i] * x[i]; 
    }
    return sqrt_gen(sum);
}


double norm2(const SpeciesJacobian& x) 
{
    double sum = double(0.0);
    for (int i  = 0; i<n_species; i++)
    {
        for (int j  = 0; j<n_species; j++)
        {
            sum += x[i][j] * x[i][j]; 
        }
    }
    return sqrt_gen(sum);
}



Species operator*(const SpeciesJacobian& A, const Species& x) 
{
    Species r = {};
    for (int i = 0; i < n_species; ++i)
    {
        for (int j = 0; j < n_species; ++j)
        {
            r[i] += A[i][j] * x[j];
        }
    }
    return r;
}


SpeciesJacobian operator*(const SpeciesJacobian& A, const SpeciesJacobian& B) 
{
    SpeciesJacobian C = {};
    for (int i = 0; i < n_species; ++i)
    {
        for (int j = 0; j < n_species; ++j)
        {
            C[i][j] = 0.0;
            for (int k = 0; k < n_species; ++k)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}


Species operator*(const double& s, const Species& a) 
{
    return scale_gen(s, a);
}


double dot(const Species& a, const Species& b) 
{
    double dot_sum = double(0.0);
    for (int i = 0; i < n_species; ++i) dot_sum += a[i]*b[i];
    return dot_sum;
}



SpeciesJacobian 
jacobian_I() 
{
    SpeciesJacobian I = {};
    for (int i = 0; i < n_species; ++i) I[i][i] = 1.0;
    return I;
}