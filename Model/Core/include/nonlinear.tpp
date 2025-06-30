#ifndef NONLINEAR_TPP
#define NONLINEAR_TPP

#include "nonlinear.h"

namespace nl
{
    // ====================================================================================================
    // Function: newtonRaphson()
    // ====================================================================================================
    template <int N,typename F1,typename F2>
    Tensor<float,N> newtonRaphson(const Tensor<float,N>& x0,F1 function,F2 jacobian,float threshold,int iterations,int subiterations)
    {
        // Initialize Variables
        int iter = 0;
        Tensor<float,N> x = x0;
        Tensor<float,N> g = function(x);
        float norm = Fastor::inner(g,g);
        float error = Fastor::max(Fastor::abs(g));
        Tensor<float,N,N> J;
        Tensor<float,N> step;
        float norm_tmp;

        // Perform Iterative Process
        while(iter < iterations && error > threshold)
        {
            // Initialize Iteration
            J = jacobian(x);
            step = Fastor::solve<SolveCompType::SimpleInv>(J,g);
            x -= step;
            g = function(x);
            norm_tmp = Fastor::inner(g,g);

            // Apply Damping
            int subiter = 0;
            while(subiter < subiterations && (norm_tmp > norm || !std::isfinite(norm_tmp)))
            {
                step /= 2;
                x += step;
                g = function(x);
                norm_tmp = Fastor::inner(g,g);
                subiter += 1;
            }

            // Conclude Iteration
            norm = norm_tmp;
            error = Fastor::max(Fastor::abs(g));
            iter += 1;
        }

        // Return
        return x;
    }

    // ====================================================================================================
    // Class: BaseNonlinear
    // ====================================================================================================
    
    // [Public] Constructor & Destructor Methods
    template <int N>
    BaseNonlinear<N>::BaseNonlinear() = default;
    template <int N>
    BaseNonlinear<N>::~BaseNonlinear() = default;
}

#endif