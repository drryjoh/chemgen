static constexpr double linear_abs_tol_default = 10 * std::numeric_limits<double>::epsilon();
static constexpr double linear_rel_tol_default = 1e-6;

static constexpr double newton_abs_tol_default = std::max(1e-8, std::numeric_limits<double>::epsilon());
static constexpr double newton_rel_tol_default = std::max(1e-8, std::numeric_limits<double>::epsilon());

static constexpr int newton_iter_max_default = 5; // backward_euler and sdirk only