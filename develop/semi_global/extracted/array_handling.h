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


ChemicalState
fill_gen_chemical_state(const double& a) 
{
    ChemicalState b;
    for (int i = 0; i<n_variables; i++)
    {
        b[i] = a;
    }
    return b;
}


Species
inv_gen(const Species& a_s) 
{
    Species b;
    for (int i = 0; i<n_species; i++)
    {
        b[i] = inv_gen(a_s[i]);
    }
    return b;
}


ChemicalState
inv_gen(ChemicalState a_s) 
{
    ChemicalState b;
    for (int i = 0; i<n_variables; i++)
    {
        b[i] = inv_gen(a_s[i]);
    }
    return b;
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

    for (int i = 0; i<n_variables; i++)
    {
        for (int j = 0; j<n_variables; j++)
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

    for (int i = 0; i<n_variables; i++)
    {
        for (int j = 0; j<n_variables; j++)
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



Species
abs(const Species& a_s) 
{
    Species b = {};
    for (int i = 0; i<n_species; i++)
    {
        b[i] = abs_gen(a_s[i]);
    }
    return b;
}


ChemicalState
abs(ChemicalState a_s) 
{
    ChemicalState b = {};
    for (int i = 0; i<n_variables; i++)
    {
        b[i] = abs_gen(a_s[i]);
    }
    return b;
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

    for (int i = 0; i<n_variables; i++)
    {
        for (int j = 0; j<n_variables; j++)
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

    for (int i = 0; i<n_variables; i++)
    {
        for (int j = 0; j<n_variables; j++)
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
set_chemical_state(const double& energy, // internal_energy or temperature
                   const Species& species) 
{
    ChemicalState y;
    for (int i = 1; i <= n_species; i++)
    {
        y[i] = species[i-1];
    }
    y[0] = energy;
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


double sum_sqr_ref(const Species& x, const Species& x_ref) 
{
    double x_ref_min = double(1);
    double x_ref_max = double(0);
    for (int i  = 0; i<n_species; i++)
    {
        if (x_ref[i] > x_ref_max)
        {
            x_ref_max = x_ref[i];
        }
        if (x_ref[i] < x_ref_min && x_ref[i] > double(1e-6))
        {
            x_ref_min = x_ref[i];
        }
    }

    double factor  = inv_gen(x_ref_max - x_ref_min);
    factor  = factor * factor;

    double sum = double(0.0);
    for (int i  = 0; i<n_species; i++)
    {
        sum += x[i] * x[i] * factor;
    }
    return sum;
}


double norm2(const ChemicalState& x) 
{
    double sum = double(0.0);
    for (int i  = 0; i<n_variables; i++)
    {
        sum += x[i] * x[i];
    }
    return sqrt_gen(sum);
}


double sum_sqr_ref(const ChemicalState& x, const ChemicalState& x_ref) 
{
    double sum = sum_sqr_ref(get_species(x), get_species(x_ref));
    sum += pow2(x[n_variables] * inv_gen(x_ref[n_variables]));
    return sum;
}


double norm2(const ChemicalState& x, const ChemicalState& x_ref) 
{
    return sqrt_gen(sum_sqr_ref(x, x_ref));
}


double norm2(const SpeciesJacobian& x) 
{
    double sum = double(0.0);
    for (int i  = 0; i < n_variables; i++)
    {
        for (int j  = 0; j < n_variables; j++)
        {
            sum += x[i][j] * x[i][j]; 
        }
    }
    return sqrt_gen(sum);
}


double error_norm(const ChemicalState& x,
                             const ChemicalState& x_ref,
                             double abs_tolerance,
                             double rel_tolerance) 
{
    return norm2(inv_gen(fill_gen_chemical_state(abs_tolerance) + scale_gen(rel_tolerance, abs(x_ref))) * x);
}



ChemicalState operator*(const SpeciesJacobian& A, const ChemicalState& x) 
{
    ChemicalState r = {};
    for (int i = 0; i < n_variables; ++i)
    {
        for (int j = 0; j < n_variables; ++j)
        {
            r[i] += A[i][j] * x[j];
        }
    }
    return r;
}


SpeciesJacobian operator*(const SpeciesJacobian& A, const SpeciesJacobian& B) 
{
    SpeciesJacobian C = {};
    for (int i = 0; i < n_variables; ++i)
    {
        for (int j = 0; j < n_variables; ++j)
        {
            C[i][j] = 0.0;
            for (int k = 0; k < n_variables; ++k)
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


ChemicalState operator*(const double& s, const ChemicalState& a) 
{
    return scale_gen(s, a);
}


double dot(const Species& a, const Species& b) 
{
    double dot_sum = double(0.0);
    for (int i = 0; i < n_species; ++i) dot_sum += a[i]*b[i];
    return dot_sum;
}


double dot(const ChemicalState& a, const ChemicalState& b) 
{
    double dot_sum = double(0.0);
    for (int i = 0; i < n_variables; ++i) dot_sum += a[i]*b[i];
    return dot_sum;
}


SpeciesJacobian 
jacobian_I() 
{
    SpeciesJacobian I = {};
    for (int i = 0; i < n_variables; ++i) I[i][i] = 1.0;
    return I;
}


ChemicalState
add_species_to_chemical_state(const ChemicalState& a_s, const Species& b_s) 
{
    ChemicalState c_s;
    c_s[0] = a_s[0];
    std::transform(a_s.begin()+1, a_s.end(), b_s.begin(), c_s.begin()+1, [](const double& a, const double& b) { return a + b; });
    return c_s;
}


Species
multiply(Species a, Species b) 
{
    return a * b;
}


Species
multiply_chain(Species a, Species a_perturbation, Species b, Species b_perturbation) 
{
    return a_perturbation * b + a * b_perturbation;
}


Species
inv(const Species& a_s) 
{
    Species b = {};
    for (int i = 0; i<n_species; i++)
    {
        b[i] = inv_gen(a_s[i]);
    }
    return b;
}


Species
operator/(const Species& a_s, const Species& b_s_divide) 
{
    Species c_s;
    Species b_s = inv(b_s_divide);
    std::transform(a_s.begin(), a_s.end(), b_s.begin(), c_s.begin(), [](const double& a, const double& b) { return a * b; });
    return c_s;
}


Species
divide(const Species& a, const Species& b) 
{
    return a / b;
}