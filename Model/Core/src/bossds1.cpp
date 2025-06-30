#include "bossds1.h"

// ====================================================================================================
// Class: BaseModel
// ====================================================================================================

// [Public] Constructor & Destructor Methods
BaseModel::BaseModel()
{
    this->fs = 44100.0f;
    this->setup();
}
BaseModel::BaseModel(float fs)
{
    this->fs = fs;
    this->setup();
}
BaseModel::~BaseModel() = default;

// [Private] Member Methods
void BaseModel::setup()
{
    this->Ts = 1/this->fs;
    this->c = 2/this->Ts;
}

// ====================================================================================================
// Class: TransistorAmplifier
// ====================================================================================================

// [Public] Constructor & Destructor Methods
TransistorAmplifier::TransistorAmplifier() : BaseModel()
{
    this->setup();
}
TransistorAmplifier::TransistorAmplifier(float fs) : BaseModel(fs)
{
    this->setup();
}
TransistorAmplifier::~TransistorAmplifier() = default;

// [Public] Other Methods
float TransistorAmplifier::processSample(float Vin)
{
    // Set Input Voltage
    this->u(0) = Vin;

    // Solve Nonlinearities
    this->p = this->Ap % this->x_1 + this->Bp % (this->u + this->u_1) + this->Cp % this->fv_1 + this->E % this->u;
    this->v = this->solveNonlinearity(this->v_1);
    this->fv = this->nonlinearFunction(this->v);

    // Update States
    this->x = this->Ax % this->x_1 + this->Bx % (this->u + this->u_1) + this->Cx % (this->fv + this->fv_1);

    // Compute Output Voltage
    this->y = this->L % this->x + this->M % this->u + this->N % this->fv;
    float Vout = this->y(0);

    // Update Memory
    this->x_1 = std::move(this->x);
    this->u_1(0) = Vin;
    this->v_1 = std::move(this->v);
    this->fv_1 = std::move(this->fv);

    // Return
    return Vout;
}

// [Private] Member Methods
void TransistorAmplifier::setup()
{
    this->Vcc = 9.0f;
    this->R6 = 100e3f;
    this->R7 = 470e3f;
    this->R8 = 10e3f;
    this->R9 = 22.0f;
    this->R10 = 100e3f;
    this->G6 = 1/this->R6;
    this->G7 = 1/this->R7;
    this->G8 = 1/this->R8;
    this->G9 = 1/this->R9;
    this->G10 = 1/this->R10;
    this->C3 = 47e-9f;
    this->C4 = 250e-12f;
    this->C5 = 68e-9f;
    this->Bf = 200.0f;
    this->Br = 0.1f;
    this->Is = 6.734e-15f;
    this->n = 1.0f;
    this->Vt = utils::THERMAL_VOLTAGE;
    this->transistor = nl::Transistor(this->Bf,this->Br,this->Is,this->n,this->Vt);
    this->x.zeros();
    this->u = {0.0f,this->Vcc};
    this->v.zeros();
    this->fv.zeros();
    this->y.zeros();
    this->p.zeros();
    this->A = {{-(this->G6 + this->G8 + this->G10)/this->C3,-(this->G8 + this->G10)/this->C3,-this->G10/this->C3},
               {-(this->G8 + this->G10)/this->C4,-(this->G7 + this->G8 + this->G10)/this->C4,-this->G10/this->C4},
               {-this->G10/this->C5,-this->G10/this->C5,-this->G10/this->C5}};
    this->B = {{(this->G6 + this->G8 + this->G10)/this->C3,-this->G8/this->C3},
               {(this->G8 + this->G10)/this->C4,-this->G8/this->C4},
               {this->G10/this->C5,0.0f}};
    this->C = {{1/this->C3,1/this->C3},
               {0.0f,1/this->C4},
               {0.0f,0.0f}};
    this->D = {{-1,0,0},
               {0,1,0}};
    this->E = {{1,0},
               {0,0}};
    this->F = {{-1/this->G9,-1/this->G9},
               {0.0f,0.0f}};
    this->L = {{-1,-1,-1}};
    this->M = {{1,0}};
    this->N = {{0,0}};
    this->Ix.eye();
    this->Iv.eye();
    this->H = Fastor::inv(this->c*this->Ix - this->A);
    this->Ax = this->H % (this->c*this->Ix + this->A);
    this->Bx = this->H % this->B;
    this->Cx = this->H % this->C;
    this->Ap = this->D % this->H % (this->c*this->Ix + this->A);
    this->Bp = this->D % this->H % this->B;
    this->Cp = this->D % this->H % this->C;
    this->K = this->D % this->H % this->C + this->F;
    this->x_dc.zeros();
    this->u_dc = {this->Vcc};
    this->v_dc.zeros();
    this->computeInitialConditions();
    this->x_1 = this->x_dc;
    this->u_1 = this->u;
    this->v_1 = this->v_dc;
    this->fv_1 = this->nonlinearFunction(this->v_1);
}
void TransistorAmplifier::computeInitialConditions()
{
    // DC System Matrices
    float Req = this->R6 + this->R7 + this->R8;
    Tensor<float,XDIM,UDIM - 1> A = {{-this->R6/Req},
                                       {-this->R7/Req},
                                       {(this->R6 + this->R7)/Req}};
    Tensor<float,XDIM,VDIM> B = {{(this->R6*this->R7 + this->R6*this->R8)/Req,this->R6*this->R8/Req},
                                   {-this->R6*this->R7/Req,this->R7*this->R8/Req},
                                   {-this->R6*this->R8/Req,-(this->R8*this->R6 + this->R8*this->R7)/Req}};
    Tensor<float,VDIM,UDIM - 1> C = {{this->R6/Req},
                                       {-this->R7/Req}};
    Tensor<float,VDIM,VDIM> D = {{-(this->R6*this->R7 + this->R6*this->R8 + this->R9*this->R6 + this->R9*this->R7 + this->R9*this->R8)/Req,-(this->R6*this->R8 + this->R9*this->R6 + this->R9*this->R7 + this->R9*this->R8)/Req},
                                   {-this->R6*this->R7/Req,this->R7*this->R8/Req}};
    
    // Define Auxiliary Functions
    auto function = [this,C,D](const Tensor<float,VDIM>& v)
    {
        return C % this->u_dc + D % this->nonlinearFunction(v) - v;
    };
    auto jacobian = [this,D](const Tensor<float,VDIM>& v)
    {
        return D % this->nonlinearJacobian(v) - this->Iv;
    };

    // Compute Initial Conditions
    this->v_dc = nl::newtonRaphson(this->v_dc,function,jacobian,this->threshold,this->iterations,this->subiterations);
    this->x_dc = A % this->u_dc + B % this->nonlinearFunction(this->v_dc);
}
Tensor<float,TransistorAmplifier::VDIM> TransistorAmplifier::nonlinearFunction(const Tensor<float,VDIM>& v)
{
    // Compute Function Value
    Tensor<float,VDIM> f = this->transistor.current(v);

    // Return
    return f;
}
Tensor<float,TransistorAmplifier::VDIM,TransistorAmplifier::VDIM> TransistorAmplifier::nonlinearJacobian(const Tensor<float,VDIM>& v)
{
    // Compute Jacobian Value
    Tensor<float,VDIM,VDIM> Jf = this->transistor.slope(v);
        
    // Return
    return Jf;
}
Tensor<float,TransistorAmplifier::VDIM> TransistorAmplifier::residualFunction(const Tensor<float,VDIM>& v)
{
    // Compute Function Value
    Tensor<float,VDIM> g = this->p + this->K % this->nonlinearFunction(v) - v;

    // Return
    return g;
}
Tensor<float,TransistorAmplifier::VDIM,TransistorAmplifier::VDIM> TransistorAmplifier::residualJacobian(const Tensor<float,VDIM>& v)
{
    // Compute Jacobian Value
    Tensor<float,VDIM,VDIM> Jg = this->K % this->nonlinearJacobian(v) - this->Iv;

    // Return
    return Jg;
}
Tensor<float,TransistorAmplifier::VDIM> TransistorAmplifier::solveNonlinearity(const Tensor<float,VDIM>& v0)
{
    // Solve Nonlinearity
    auto function = [this](const Tensor<float,VDIM>& v)
    {
        return this->residualFunction(v);
    };
    auto jacobian = [this](const Tensor<float,VDIM>& v)
    {
        return this->residualJacobian(v);
    };
    Tensor<float,VDIM> v = nl::newtonRaphson(v0,function,jacobian,this->threshold,this->iterations,this->subiterations);

    // Return
    return v;
}

// ====================================================================================================
// Class: OperationalAmplifier
// ====================================================================================================

// [Public] Constructor & Destructor Methods
OperationalAmplifier::OperationalAmplifier() : BaseModel()
{
    this->distortion = 0.5f;
    this->setup();
}
OperationalAmplifier::OperationalAmplifier(float fs,float distortion) : BaseModel(fs)
{
    this->distortion = utils::limit(distortion,0,1,utils::POT_MARGIN);
    this->setup();
}
OperationalAmplifier::~OperationalAmplifier() = default;

// [Public] Other Methods
float OperationalAmplifier::processSample(float Vin)
{
    // Compute Feedback Current
    this->Vin_port.setVoltage(Vin);
    this->series_root.collectIncidentWaves();
    this->series_root.propagateReflectedWaves();
    float Ifb = -this->Vin_port.waveToCurrent();

    // Compute Feedback Voltage
    this->Ifb_port.setCurrent(Ifb);
    this->parallel_root.collectIncidentWaves();
    this->parallel_root.propagateReflectedWaves();
    float Vfb = this->Ifb_port.waveToVoltage();

    // Compute Output Voltage
    float Vout = utils::limit(Vfb + Vin,-this->Vsat,this->Vsat);

    // Return
    return Vout;
}
void OperationalAmplifier::setDistortion(float distortion)
{
    // Check Change
    if(this->distortion != distortion)
    {
        // Update Distortion Value
        this->distortion = utils::limit(distortion,0,1,utils::POT_MARGIN);

        // Update WDF Structure
        this->updateStructure();
    }
}
float OperationalAmplifier::getDistortion()
{
    return this->distortion;
}

// [Private] Member Methods
void OperationalAmplifier::setup()
{
    this->R11 = 100e3f;
    this->R13 = 4.7e3f;
    this->VR1 = 100e3f;
    this->C7 = 100e-12f;
    this->C8 = 0.47e-6f;
    this->Vsat = 4.5f;
    this->Ra = this->VR1*(1 - this->distortion);
    this->Rb = this->VR1*this->distortion;
    this->Vin_port = wdf::ResistiveVoltageSource(this->Ra);
    this->R13_port = wdf::Resistor(this->R13);
    this->C8_port = wdf::Capacitor(this->C8,this->fs);
    this->series_root = wdf::SeriesAdaptorRoot<wdf::ResistiveVoltageSource,wdf::Resistor,wdf::Capacitor>(this->Vin_port,this->R13_port,this->C8_port);
    this->Ifb_port = wdf::ResistiveCurrentSource(this->Rb);
    this->C7_port = wdf::Capacitor(this->C7,this->fs);
    this->parallel_root = wdf::ParallelAdaptorRoot<wdf::ResistiveCurrentSource,wdf::Capacitor>(this->Ifb_port,this->C7_port);
}
void OperationalAmplifier::updateStructure()
{
    // Update Variable Resistors
    this->Ra = this->VR1*(1 - this->distortion);
    this->Rb = this->VR1*this->distortion;

    // Update Port Resistances
    this->Vin_port.setPortResistance(this->Ra);
    this->Ifb_port.setPortResistance(this->Rb);

    // Update Scattering Matrices
    this->series_root.computeScatteringMatrix();
    this->parallel_root.computeScatteringMatrix();
}

// ====================================================================================================
// Class: DiodeClipper
// ====================================================================================================

// [Public] Constructor & Destructor Methods
DiodeClipper::DiodeClipper() : BaseModel()
{
    this->setup();
}
DiodeClipper::DiodeClipper(float fs) : BaseModel(fs)
{
    this->setup();
}
DiodeClipper::~DiodeClipper() = default;

// [Public] Other Methods
float DiodeClipper::processSample(float Vin)
{
    // Set Input Voltage
    this->u(0) = Vin;

    // Solve Nonlinearities
    this->p = this->Aw % this->x + this->Cw % this->u;
    this->w = this->solveNonlinearity(this->w);
    this->z = this->nonlinearFunction(this->w);

    // Compute Output Voltage
    this->y = this->Ay % this->x + this->By % this->z + this->Cy % this->u;
    float Vout = -this->y(last);

    // Update States
    this->dx = this->Ax % this->x + this->Bx % this->z + this->Cx % this->u;
    this->x += this->dx;
    
    // Return
    return Vout;
}

// [Private] Member Methods
void DiodeClipper::setup()
{
    this->R14 = 2.2e3f;
    this->G14 = 1/this->R14;
    this->C9 = 0.47e-6f;
    this->C10 = 0.01e-6f;
    this->Is = 2.52e-9f;
    this->n = 1.752f;
    this->diodes = nl::Diodes(this->Is,this->n,this->Vt);
    this->x.zeros();
    this->dx.zeros();
    this->w.zeros();
    this->z.zeros();
    this->u.zeros();
    this->y.zeros();
    this->p.zeros();
    this->permutation = {3,4,2,0,5,1,6};
    this->incident_matrix = {{0,-1,0,0,0,-1,-1},
                             {0,0,1,0,0,1,0},
                             {0,0,-1,1,0,0,0},
                             {1,1,0,-1,1,0,0},
                             {-1,0,0,0,-1,0,1}};
    this->computeMatrixJ();
    this->Jx = this->J(fseq<0,SDIM>(),fseq<0,SDIM>());
    this->Jw = this->J(fseq<SDIM,SDIM + DDIM>(),fseq<SDIM,SDIM + DDIM>());
    this->Jy = this->J(fseq<SDIM + DDIM,SDIM + DDIM + PDIM>(),fseq<SDIM + DDIM,SDIM + DDIM + PDIM>());
    this->K = -this->J(fseq<0,SDIM>(),fseq<SDIM,SDIM + DDIM>());
    this->Gx = -this->J(fseq<0,SDIM>(),fseq<SDIM + DDIM,SDIM + DDIM + PDIM>());
    this->Gw = -this->J(fseq<SDIM,SDIM + DDIM>(),fseq<SDIM + DDIM,SDIM + DDIM + PDIM>());
    this->Ix.eye();
    this->Iw.eye();
    this->Q.zeros();
    Fastor::diag(this->Q) = Tensor<float,SDIM>{1/this->C9,1/this->C10};
    this->D = Fastor::inv(this->Ix/this->Ts - (this->Jx % this->Q)/2);
    this->Ax = this->D % this->Jx % this->Q;
    this->Bx = -this->D % this->K;
    this->Cx = -this->D % this->Gx;
    this->Aw = (Fastor::transpose(this->K) % this->Q % (2*this->Ix + this->Ax))/2;
    this->Bw = this->Jw + (Fastor::transpose(this->K) % this->Q % this->Bx)/2;
    this->Cw = -this->Gw + (Fastor::transpose(this->K) % this->Q % this->Cx)/2;
    this->Ay = -((Fastor::transpose(this->Gx) % this->Q % (2*this->Ix + this->Ax))/2);
    this->By = -(Fastor::transpose(this->Gw) + (Fastor::transpose(this->Gx) % this->Q % this->Bx)/2);
    this->Cy = -(this->Jy + (Fastor::transpose(this->Gx) % this->Q % this->Cx)/2);
}
void DiodeClipper::computeMatrixJ()
{
    // Compute Auxiliary Matrices
    Tensor<float,BDIM,BDIM> P;
    P.eye();
    P = P(this->permutation,all);
    Tensor<float,NDIM - 1,BDIM - NDIM + 1> g1 = this->incident_matrix(fseq<1,last>(),fseq<0,BDIM - NDIM + 1>());
    Tensor<float,NDIM - 1,NDIM - 1> g2 = this->incident_matrix(fseq<1,last>(),fseq<BDIM - NDIM + 1,last>());
    Tensor<float,NDIM - 1,BDIM - NDIM + 1> g = Fastor::inv(g2) % g1;

    // Compute Matrix J
    Tensor<float,BDIM,BDIM> J_tmp;
    J_tmp.zeros();
    J_tmp(seq(BDIM - g.dimension(0),last),seq(first,g.dimension(1))) = -g;
    J_tmp(seq(first,g.dimension(1)),seq(BDIM - g.dimension(0),last)) = Fastor::transpose(g);
    J_tmp = P % J_tmp % Fastor::transpose(P);
    this->J = J_tmp(fseq<first,last - 1>(),fseq<first,last - 1>());
}
Tensor<float,DiodeClipper::DDIM> DiodeClipper::nonlinearFunction(const Tensor<float,DDIM>& w)
{
    // Compute Function Value
    Tensor<float,DDIM> z = {w(0)/this->R14,this->diodes.current(w(fseq<1,last>()))(0)};

    // Return
    return z;
}
Tensor<float,DiodeClipper::DDIM,DiodeClipper::DDIM> DiodeClipper::nonlinearJacobian(const Tensor<float,DDIM>& w)
{
    // Compute Jacobian Value
    Tensor<float,DDIM,DDIM> Jz = {{1/this->R14,0.0f},{0.0f,this->diodes.slope(w(fseq<1,last>()))(0,0)}};
        
    // Return
    return Jz;
}
Tensor<float,DiodeClipper::DDIM> DiodeClipper::residualFunction(const Tensor<float,DDIM>& w)
{
    // Compute Function Value
    Tensor<float,DDIM> g = this->p + this->Bw % this->nonlinearFunction(w) - w;

    // Return
    return g;
}
Tensor<float,DiodeClipper::DDIM,DiodeClipper::DDIM> DiodeClipper::residualJacobian(const Tensor<float,DDIM>& w)
{
    // Compute Jacobian Value
    Tensor<float,DDIM,DDIM> Jg = this->Bw % this->nonlinearJacobian(w) - this->Iw;

    // Return
    return Jg;
}
Tensor<float,DiodeClipper::DDIM> DiodeClipper::solveNonlinearity(const Tensor<float,DDIM>& w0)
{
    // Solve Nonlinearity
    auto function = [this](const Tensor<float,DDIM>& w)
    {
        return this->residualFunction(w);
    };
    auto jacobian = [this](const Tensor<float,DDIM>& w)
    {
        return this->residualJacobian(w);
    };
    Tensor<float,DDIM> w = nl::newtonRaphson(w0,function,jacobian,this->threshold,this->iterations,this->subiterations);

    // Return
    return w;
}

// ====================================================================================================
// Class: PassiveFilter
// ====================================================================================================

// [Public] Constructor & Destructor Methods
PassiveFilter::PassiveFilter() : BaseModel()
{
    this->tone = 0.5f;
    this->setup();
}
PassiveFilter::PassiveFilter(float fs,float tone) : BaseModel(fs)
{
    this->tone = tone;
    this->setup();
}
PassiveFilter::~PassiveFilter() = default;

// [Public] Other Methods
float PassiveFilter::processSample(float Vin)
{
    // Set Input Voltage
    this->x(0) = Vin;

    // Compute Output Voltage
    float Vout = Fastor::inner(this->b,this->x) - Fastor::inner(this->a,this->y);

    // Update Memory
    for(int n = ORDER; n > 0; n--)
    {
        this->x(n) = this->x(n - 1);
        if(n > 1)
        {
            this->y(n - 1) = this->y(n - 2);
        }
    }
    y(0) = Vout;

    // Return
    return Vout;
}
void PassiveFilter::setTone(float tone)
{
    // Check Change
    if(this->tone != tone)
    {
        // Update Tone Value
        this->tone = utils::limit(tone,0,1);

        // Update Coefficients
        this->updateCoeffs();
    }
}
float PassiveFilter::getTone()
{
    return this->tone;
}

// [Private] Member Methods
void PassiveFilter::setup()
{
    this->R15 = 2.2e3f;
    this->R16 = 6.8e3f;
    this->R17 = 6.8e3f;
    this->VR2 = 20e3f;
    this->VR3 = 100e3f;
    this->C11 = 22e-9f;
    this->C12 = 0.1e-6f;
    this->x.zeros();
    this->y.zeros();
    this->exponents.arange(0);
    this->initializeCoeffs();
}
void PassiveFilter::initializeCoeffs()
{
    // Coefficient B0
    float k_B0_0 = -this->C11*this->R15*this->R17*this->VR3*this->c - this->C11*this->R15*this->VR3*this->VR2*this->c - this->C11*this->R16*this->R17*this->VR3*this->c - this->C11*this->R17*this->VR3*this->VR2*this->c - this->R17*this->VR3 - this->VR3*this->VR2;
    float k_B0_1 = -this->C11*this->C12*this->R16*this->R17*this->VR3*this->VR2*std::pow(this->c,2) + this->C11*this->R15*this->VR3*this->VR2*this->c + this->VR3*this->VR2;
    float k_B0_2 = 0.0f;
    this->k_B0 = {k_B0_0,k_B0_1,k_B0_2};

    // Coefficient B1
    float k_B1_0 = -2*this->R17*this->VR3 - 2*this->VR3*this->VR2;
    float k_B1_1 = 2*this->C11*this->C12*this->R16*this->R17*this->VR3*this->VR2*std::pow(this->c,2) + 2*this->VR3*this->VR2;
    float k_B1_2 = 0.0f;
    this->k_B1 = {k_B1_0,k_B1_1,k_B1_2};

    // Coefficient B2
    float k_B2_0 = this->C11*this->R15*this->R17*this->VR3*this->c + this->C11*this->R15*this->VR3*this->VR2*this->c + this->C11*this->R16*this->R17*this->VR3*this->c + this->C11*this->R17*this->VR3*this->VR2*this->c - this->R17*this->VR3 - this->VR3*this->VR2;
    float k_B2_1 = -this->C11*this->C12*this->R16*this->R17*this->VR3*this->VR2*std::pow(this->c,2) - this->C11*this->R15*this->VR3*this->VR2*this->c + this->VR3*this->VR2;
    float k_B2_2 = 0.0f;
    this->k_B2 = {k_B2_0,k_B2_1,k_B2_2};

    // Coefficient Α0
    float k_A0_0 = -this->C11*this->C12*this->R15*this->R16*this->R17*this->VR3*std::pow(this->c,2) - this->C11*this->C12*this->R15*this->R16*this->VR3*this->VR2*std::pow(this->c,2) - this->C11*this->C12*this->R16*this->R17*this->VR3*this->VR2*std::pow(this->c,2) - this->C11*this->R15*this->R16*this->R17*this->c - this->C11*this->R15*this->R16*this->VR3*this->c - this->C11*this->R15*this->R16*this->VR2*this->c - this->C11*this->R15*this->R17*this->VR3*this->c - this->C11*this->R15*this->VR3*this->VR2*this->c - this->C11*this->R16*this->R17*this->VR3*this->c - this->C11*this->R16*this->R17*this->VR2*this->c - this->C11*this->R17*this->VR3*this->VR2*this->c - this->C12*this->R16*this->R17*this->VR3*this->c - this->C12*this->R16*this->VR3*this->VR2*this->c - this->R16*this->R17 - this->R16*this->VR3 - this->R16*this->VR2 - this->R17*this->VR3 - this->VR3*this->VR2;
    float k_A0_1 = -this->C11*this->C12*this->R15*this->R16*this->R17*this->VR2*std::pow(this->c,2) - this->C11*this->C12*this->R15*this->R16*std::pow(this->VR2,2)*std::pow(this->c,2) - this->C11*this->C12*this->R16*this->R17*std::pow(this->VR2,2)*std::pow(this->c,2) + this->C11*this->R15*this->R16*this->VR2*this->c - this->C11*this->R15*this->R17*this->VR2*this->c - this->C11*this->R15*std::pow(this->VR2,2)*this->c + this->C11*this->R16*this->R17*this->VR2*this->c - this->C11*this->R17*std::pow(this->VR2,2)*this->c - this->C12*this->R16*this->R17*this->VR2*this->c - this->C12*this->R16*std::pow(this->VR2,2)*this->c + this->R16*this->VR2 - this->R17*this->VR2 - std::pow(this->VR2,2);
    float k_A0_2 = this->C11*this->C12*this->R15*this->R16*std::pow(this->VR2,2)*std::pow(this->c,2) + this->C11*this->C12*this->R16*this->R17*std::pow(this->VR2,2)*std::pow(this->c,2) + this->C11*this->R15*std::pow(this->VR2,2)*this->c + this->C11*this->R17*std::pow(this->VR2,2)*this->c + this->C12*this->R16*std::pow(this->VR2,2)*this->c + std::pow(this->VR2,2);
    this->k_A0 = {k_A0_0,k_A0_1,k_A0_2};
    
    // Coefficient Α1
    float k_A1_0 = 2*this->C11*this->C12*this->R15*this->R16*this->R17*this->VR3*std::pow(this->c,2) + 2*this->C11*this->C12*this->R15*this->R16*this->VR3*this->VR2*std::pow(this->c,2) + 2*this->C11*this->C12*this->R16*this->R17*this->VR3*this->VR2*std::pow(this->c,2) - 2*this->R16*this->R17 - 2*this->R16*this->VR3 - 2*this->R16*this->VR2 - 2*this->R17*this->VR3 - 2*this->VR3*this->VR2;
    float k_A1_1 = 2*this->C11*this->C12*this->R15*this->R16*this->R17*this->VR2*std::pow(this->c,2) + 2*this->C11*this->C12*this->R15*this->R16*std::pow(this->VR2,2)*std::pow(this->c,2) + 2*this->C11*this->C12*this->R16*this->R17*std::pow(this->VR2,2)*std::pow(this->c,2) + 2*this->R16*this->VR2 - 2*this->R17*this->VR2 - 2*std::pow(this->VR2,2);
    float k_A1_2 = -2*this->C11*this->C12*this->R15*this->R16*std::pow(this->VR2,2)*std::pow(this->c,2) - 2*this->C11*this->C12*this->R16*this->R17*std::pow(this->VR2,2)*std::pow(this->c,2) + 2*std::pow(this->VR2,2);
    this->k_A1 = {k_A1_0,k_A1_1,k_A1_2};
    
    // Coefficient Α2
    float k_A2_0 = -this->C11*this->C12*this->R15*this->R16*this->R17*this->VR3*std::pow(this->c,2) - this->C11*this->C12*this->R15*this->R16*this->VR3*this->VR2*std::pow(this->c,2) - this->C11*this->C12*this->R16*this->R17*this->VR3*this->VR2*std::pow(this->c,2) + this->C11*this->R15*this->R16*this->R17*this->c + this->C11*this->R15*this->R16*this->VR3*this->c + this->C11*this->R15*this->R16*this->VR2*this->c + this->C11*this->R15*this->R17*this->VR3*this->c + this->C11*this->R15*this->VR3*this->VR2*this->c + this->C11*this->R16*this->R17*this->VR3*this->c + this->C11*this->R16*this->R17*this->VR2*this->c + this->C11*this->R17*this->VR3*this->VR2*this->c + this->C12*this->R16*this->R17*this->VR3*this->c + this->C12*this->R16*this->VR3*this->VR2*this->c - this->R16*this->R17 - this->R16*this->VR3 - this->R16*this->VR2 - this->R17*this->VR3 - this->VR3*this->VR2;
    float k_A2_1 = -this->C11*this->C12*this->R15*this->R16*this->R17*this->VR2*std::pow(this->c,2) - this->C11*this->C12*this->R15*this->R16*std::pow(this->VR2,2)*std::pow(this->c,2) - this->C11*this->C12*this->R16*this->R17*std::pow(this->VR2,2)*std::pow(this->c,2) - this->C11*this->R15*this->R16*this->VR2*this->c + this->C11*this->R15*this->R17*this->VR2*this->c + this->C11*this->R15*std::pow(this->VR2,2)*this->c - this->C11*this->R16*this->R17*this->VR2*this->c + this->C11*this->R17*std::pow(this->VR2,2)*this->c + this->C12*this->R16*this->R17*this->VR2*this->c + this->C12*this->R16*std::pow(this->VR2,2)*this->c + this->R16*this->VR2 - this->R17*this->VR2 - std::pow(this->VR2,2);
    float k_A2_2 = this->C11*this->C12*this->R15*this->R16*std::pow(this->VR2,2)*std::pow(this->c,2) + this->C11*this->C12*this->R16*this->R17*std::pow(this->VR2,2)*std::pow(this->c,2) - this->C11*this->R15*std::pow(this->VR2,2)*this->c - this->C11*this->R17*std::pow(this->VR2,2)*this->c - this->C12*this->R16*std::pow(this->VR2,2)*this->c + std::pow(this->VR2,2);
    this->k_A2 = {k_A2_0,k_A2_1,k_A2_2};

    // Coefficient Initialization
    this->updateCoeffs();
}
void PassiveFilter::updateCoeffs()
{
    // Coefficient Evaluation
    Fastor::Tensor<float,ORDER + 1> powers;
    powers.fill(this->tone);
    powers = Fastor::pow(powers,this->exponents);
    float B0 = Fastor::inner(this->k_B0,powers);
    float B1 = Fastor::inner(this->k_B1,powers);
    float B2 = Fastor::inner(this->k_B2,powers);
    float A0 = Fastor::inner(this->k_A0,powers);
    float A1 = Fastor::inner(this->k_A1,powers);
    float A2 = Fastor::inner(this->k_A2,powers);

    // Coefficient Normalization
    this->b = {B0,B1,B2};
    this->b = this->b/A0;
    this->a = {A1,A2};
    this->a = this->a/A0;
}

// ====================================================================================================
// Class: BossDS1
// ====================================================================================================

// [Public] Constructor & Destructor Methods
BossDS1::BossDS1() : BaseModel(),
                     transistor_amplifier(fs),
                     operational_amplifier(fs),
                     diode_clipper(fs),
                     passive_filter(fs)
{
    this->distortion = this->operational_amplifier.getDistortion();
    this->tone = this->passive_filter.getTone();
    this->level = 1.0f;
}
BossDS1::BossDS1(float fs,float distortion,float tone,float level) : BaseModel(fs),
                                                                     transistor_amplifier(fs),
                                                                     operational_amplifier(fs,distortion),
                                                                     diode_clipper(fs),
                                                                     passive_filter(fs,tone)
{
    this->distortion = this->operational_amplifier.getDistortion();
    this->tone = this->passive_filter.getTone();
    this->level = utils::limit(level,0,1);
}
BossDS1::~BossDS1() = default;

// [Public] Other Methods
float BossDS1::processSample(float Vin)
{
    // Compute Output Voltage
    float Vout = Vin;
    Vout = this->transistor_amplifier.processSample(Vout);
    Vout = this->operational_amplifier.processSample(Vout);
    Vout = this->diode_clipper.processSample(Vout);
    Vout = this->passive_filter.processSample(Vout);
    Vout = this->level*Vout;

    // Return
    return Vout;
}
void BossDS1::setDistortion(float distortion)
{
    // Update Distortion Value
    this->operational_amplifier.setDistortion(distortion);
}
void BossDS1::setTone(float tone)
{
    // Update Tone Value
    this->passive_filter.setTone(tone);
}
void BossDS1::setLevel(float level)
{
    // Check Change
    if(this->level != level)
    {
        // Update Level Value
        this->level = utils::limit(level,0,1);
    }
}
float BossDS1::getDistortion()
{
    return this->distortion;
}
float BossDS1::getTone()
{
    return this->tone;
}
float BossDS1::getLevel()
{
    return this->level;
}