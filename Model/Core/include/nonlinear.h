#ifndef NONLINEAR_H
#define NONLINEAR_H

#include <cmath>

#include "utilities.h"

#include <Fastor/Fastor.h>

using namespace Fastor;

namespace nl
{
    // ====================================================================================================
    // Function: newtonRaphson()
    // ====================================================================================================
    template <int N,typename F1,typename F2>
    Tensor<float,N> newtonRaphson(const Tensor<float,N>& x0,F1 function,F2 jacobian,float threshold = utils::THRESHOLD,int iterations = utils::ITERATIONS,int subiterations = utils::SUBITERATIONS);

    // ====================================================================================================
    // Class: BaseNonlinear
    // ====================================================================================================
    template <int N>
    class BaseNonlinear
    {
        public:
        // Constructor & Destructor Methods
        BaseNonlinear();
        virtual ~BaseNonlinear();
        
        // Other Methods
        virtual Tensor<float,N> current(const Tensor<float,N>& v) = 0;
        virtual Tensor<float,N,N> slope(const Tensor<float,N>& v) = 0;
        
        protected:
        // Static Variables
        static constexpr int DIM = N;
    };

    // ====================================================================================================
    // Class: Diodes
    // ====================================================================================================
    class Diodes final : public BaseNonlinear<1>
    {
        public:
        // Constructor & Destructor Methods
        Diodes();
        Diodes(float Is,float n = 1,float Vt = utils::THERMAL_VOLTAGE);
        ~Diodes();
    
        // Other Methods
        Tensor<float,Diodes::DIM> current(const Tensor<float,Diodes::DIM>& v) override;
        Tensor<float,Diodes::DIM,Diodes::DIM> slope(const Tensor<float,Diodes::DIM>& v) override;

        private:
        // Member Variables
        float Is;
        float n;
        float Vt;
        float k;
        float d;
        float kd;

        // Member Methods
        void setup();
    };

    // ====================================================================================================
    // Class: Transistor
    // ====================================================================================================
    class Transistor final : public BaseNonlinear<2>
    {
        public:
        // Constructor & Destructor Methods
        Transistor();
        Transistor(float Bf,float Br,float Is,float n = 1,float Vt = utils::THERMAL_VOLTAGE);
        ~Transistor();
    
        // Other Methods
        Tensor<float,Transistor::DIM> current(const Tensor<float,Transistor::DIM>& v) override;
        Tensor<float,Transistor::DIM,Transistor::DIM> slope(const Tensor<float,Transistor::DIM>& v) override;

        private:
        // Member Variables
        float Bf;
        float Br;
        float Is;
        float n;
        float Vt;
        float k;
        float b0;
        float b1;
        float c0;
        float c1;
        float kb0;
        float kb1;
        float kc0;
        float kc1;

        // Member Methods
        void setup();
    };
}

#include "nonlinear.tpp"

#endif