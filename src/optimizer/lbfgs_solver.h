#pragma once
#include "types.h"
#include <functional>

struct LBFGSConfig {
    int max_iter = 200; // reduced: good init+analytic grad converges faster
    int history_size = 10;
    int max_ls_iter = 20;
    double grad_tol = 1e-4; // loosened: curve quality doesn't need 1e-5
    double func_tol = 1e-7;
    double wolfe_c1 = 1e-4;
    double wolfe_c2 = 0.9;
};

struct SolveResult {
    VecXd x;
    double final_cost = 0;
    int iterations = 0;
    bool converged = false;
    int fn_evals = 0; // diagnostic: total cost-function evaluations
};

// Cost function: returns f, fills grad.
// The function MUST compute grad correctly — lineSearch no longer re-evaluates.
using CostFn = std::function<double(const VecXd&, VecXd&)>;

class LBFGSSolver {
public:
    explicit LBFGSSolver(const LBFGSConfig& cfg = {}) : cfg_(cfg) {}

    SolveResult solve(CostFn fn, const VecXd& x0);
    SolveResult solveWarm(CostFn fn, const VecXd& x_warm);

private:
    LBFGSConfig cfg_;

    // Strong Wolfe line search.
    // Returns accepted alpha and stores accepted (xn, fn_val, gn).
    // Does NOT recompute fn at the accepted point — reuses the last evaluation.
    double lineSearch(CostFn& fn, const VecXd& x, const VecXd& dir,
                      double f0, const VecXd& g0, VecXd& xn, double& fn_val, VecXd& gn, int& evals);
};
