#include "nonlinear.h"

namespace nl
{
    // ====================================================================================================
    // Class: Diodes
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    Diodes::Diodes()
    {
        this->Is = 1e-9f;
        this->n = 1.0f;
        this->Vt = utils::THERMAL_VOLTAGE;
        this->setup();
    }
    Diodes::Diodes(float Is,float n,float Vt)
    {
        this->Is = Is;
        this->n = n;
        this->Vt = Vt;
        this->setup();
    }
    Diodes::~Diodes() = default;

    // [Public] Other Methods
    Tensor<float,Diodes::DIM> Diodes::current(const Tensor<float,Diodes::DIM>& v)
    {
        // Calculate Current
        float i = this->d*std::sinh(v(0)*this->k);
        Tensor<float,Diodes::DIM> current = {i};

        // Return
        return current;
    }
    Tensor<float,Diodes::DIM,Diodes::DIM> Diodes::slope(const Tensor<float,Diodes::DIM>& v)
    {
        // Calculate Slope
        float s = this->kd*std::cosh(v(0)*this->k);
        Tensor<float,Diodes::DIM,Diodes::DIM> slope = {{s}};

        // Return
        return slope;
    }

    // [Private] Member Methods
    void Diodes::setup()
    {
        this->k = 1/(this->n*this->Vt);
        this->d = 2*this->Is;
        this->kd = this->k*this->d;
    }

    // ====================================================================================================
    // Class: Transistor
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    Transistor::Transistor()
    {
        this->Bf = 100.0f;
        this->Br = 0.1f;
        this->Is = 1e-9f;
        this->n = 1.0f;
        this->Vt = utils::THERMAL_VOLTAGE;
        this->setup();
    }
    Transistor::Transistor(float Bf,float Br,float Is,float n,float Vt)
    {
        this->Bf = Bf;
        this->Br = Br;
        this->Is = Is;
        this->n = n;
        this->Vt = Vt;
        this->setup();
    }
    Transistor::~Transistor() = default;

    // [Public] Other Methods
    Tensor<float,Transistor::DIM> Transistor::current(const Tensor<float,Transistor::DIM>& v)
    {
        // Calculate Current
        float exp0 = std::exp(v(0)*this->k) - 1;
        float exp1 = std::exp(v(1)*this->k) - 1;
        float ib = this->b0*exp0 + this->b1*exp1;
        float ic = this->c0*exp0 + this->c1*exp1;
        Tensor<float,Transistor::DIM> current = {ib,ic};

        // Return
        return current;
    }
    Tensor<float,Transistor::DIM,Transistor::DIM> Transistor::slope(const Tensor<float,Transistor::DIM>& v)
    {
        // Calculate Slope
        float exp0 = std::exp(v(0)*this->k);
        float exp1 = std::exp(v(1)*this->k);
        float sb0 = this->kb0*exp0;
        float sb1 = this->kb1*exp1;
        float sc0 = this->kc0*exp0;
        float sc1 = this->kc1*exp1;
        Tensor<float,Transistor::DIM,Transistor::DIM> slope = {{sb0,sb1},
                                                               {sc0,sc1}};

        // Return
        return slope;
    }

    // [Private] Member Methods
    void Transistor::setup()
    {
        this->k = 1/(this->n*this->Vt);
        this->b0 = this->Is/this->Bf;
        this->b1 = this->Is/this->Br;
        this->c0 = this->Is;
        this->c1 = -this->Is/(this->Br/(1 + this->Br));
        this->kb0 = this->k*this->b0;
        this->kb1 = this->k*this->b1;
        this->kc0 = this->k*this->c0;
        this->kc1 = this->k*this->c1;
    }
}