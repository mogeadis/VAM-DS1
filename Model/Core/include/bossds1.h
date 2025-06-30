#ifndef BOSSDS1_H
#define BOSSDS1_H

#include <cmath>

#include "wdf.h"
#include "nonlinear.h"
#include "utilities.h"

#include <Fastor/Fastor.h>

using namespace Fastor;

// ====================================================================================================
// Class: BaseModel
// ====================================================================================================
class BaseModel
{
    public:
    // Constructor & Destructor Methods
    BaseModel();
    BaseModel(float fs);
    virtual ~BaseModel();

    // Other Methods
    virtual float processSample(float Vin) = 0;

    protected:
    // Member Variables
    float fs;
    float Ts;
    float c;

    private:
    // Member Methods
    void setup();
};

// ====================================================================================================
// Class: TransistorAmplifier
// ====================================================================================================
class TransistorAmplifier final : public BaseModel
{
    public:
    // Constructor & Destructor Methods
    TransistorAmplifier();
    TransistorAmplifier(float fs);
    ~TransistorAmplifier();

    // Other Methods
    float processSample(float Vin) override;

    private:
    // Static Variables
    static constexpr int XDIM = 3;
    static constexpr int UDIM = 2;
    static constexpr int VDIM = 2;
    static constexpr int YDIM = 1;

    // Member Variables
    float Vcc;
    float R6;
    float R7;
    float R8;
    float R9;
    float R10;
    float G6;
    float G7;
    float G8;
    float G9;
    float G10;
    float C3;
    float C4;
    float C5;
    float Bf;
    float Br;
    float Is;
    float n;
    float Vt = utils::THERMAL_VOLTAGE;
    nl::Transistor transistor;
    float threshold = utils::THRESHOLD;
    int iterations = utils::ITERATIONS;
    int subiterations = utils::SUBITERATIONS;
    Tensor<float,XDIM> x;
    Tensor<float,UDIM> u;
    Tensor<float,VDIM> v;
    Tensor<float,VDIM> fv;
    Tensor<float,YDIM> y;
    Tensor<float,VDIM> p;
    Tensor<float,XDIM,XDIM> A;
    Tensor<float,XDIM,UDIM> B;
    Tensor<float,XDIM,VDIM> C;
    Tensor<float,VDIM,XDIM> D;
    Tensor<float,VDIM,UDIM> E;
    Tensor<float,VDIM,VDIM> F;
    Tensor<float,YDIM,XDIM> L;
    Tensor<float,YDIM,UDIM> M;
    Tensor<float,YDIM,VDIM> N;
    Tensor<float,XDIM,XDIM> Ix;
    Tensor<float,VDIM,VDIM> Iv;
    Tensor<float,XDIM,XDIM> H;
    Tensor<float,XDIM,XDIM> Ax;
    Tensor<float,XDIM,UDIM> Bx;
    Tensor<float,XDIM,VDIM> Cx;
    Tensor<float,VDIM,XDIM> Ap;
    Tensor<float,VDIM,UDIM> Bp;
    Tensor<float,VDIM,VDIM> Cp;
    Tensor<float,VDIM,VDIM> K;
    Tensor<float,XDIM> x_dc;
    Tensor<float,UDIM - 1> u_dc;
    Tensor<float,VDIM> v_dc;
    Tensor<float,XDIM> x_1;
    Tensor<float,UDIM> u_1;
    Tensor<float,VDIM> v_1;
    Tensor<float,VDIM> fv_1;

    // Member Methods
    void setup();
    void computeInitialConditions();
    Tensor<float,VDIM> nonlinearFunction(const Tensor<float,VDIM>& w);
    Tensor<float,VDIM,VDIM> nonlinearJacobian(const Tensor<float,VDIM>& w);
    Tensor<float,VDIM> residualFunction(const Tensor<float,VDIM>& w);
    Tensor<float,VDIM,VDIM> residualJacobian(const Tensor<float,VDIM>& w);
    Tensor<float,VDIM> solveNonlinearity(const Tensor<float,VDIM>& w0);
};

// ====================================================================================================
// Class: OperationalAmplifier
// ====================================================================================================
class OperationalAmplifier final : public BaseModel
{
    public:
    // Constructor & Destructor Methods
    OperationalAmplifier();
    OperationalAmplifier(float fs,float distortion = 0.5f);
    ~OperationalAmplifier();

    // Other Methods
    float processSample(float Vin) override;
    void setDistortion(float distortion);
    float getDistortion();

    private:
    // Member Variables
    float distortion;
    float R11;
    float R13;
    float VR1;
    float C7;
    float C8;
    float Vsat;
    float Ra;
    float Rb;
    wdf::ResistiveVoltageSource Vin_port;
    wdf::Resistor R13_port;
    wdf::Capacitor C8_port;
    wdf::SeriesAdaptorRoot<wdf::ResistiveVoltageSource,wdf::Resistor,wdf::Capacitor> series_root;
    wdf::ResistiveCurrentSource Ifb_port;
    wdf::Capacitor C7_port;
    wdf::ParallelAdaptorRoot<wdf::ResistiveCurrentSource,wdf::Capacitor> parallel_root;

    // Member Methods
    void setup();
    void updateStructure();
};

// ====================================================================================================
// Class: DiodeClipper
// ====================================================================================================
class DiodeClipper final : public BaseModel
{
    public:
    // Constructor & Destructor Methods
    DiodeClipper();
    DiodeClipper(float fs);
    ~DiodeClipper();

    // Other Methods
    float processSample(float Vin) override;

    private:
    // Static Variables
    static constexpr int SDIM = 2;
    static constexpr int DDIM = 2;
    static constexpr int PDIM = 3 - 1;
    static constexpr int NDIM = 4 + 1;
    static constexpr int BDIM = SDIM + DDIM + PDIM + 1;

    // Member Variables
    float R14;
    float G14;
    float C9;
    float C10;
    float Is;
    float n;
    float Vt = utils::THERMAL_VOLTAGE;
    nl::Diodes diodes;
    float threshold = utils::THRESHOLD;
    int iterations = utils::ITERATIONS;
    int subiterations = utils::SUBITERATIONS;
    Tensor<float,SDIM> x;
    Tensor<float,SDIM> dx;
    Tensor<float,DDIM> w;
    Tensor<float,DDIM> z;
    Tensor<float,PDIM> u;
    Tensor<float,PDIM> y;
    Tensor<float,DDIM> p;
    Tensor<int,BDIM> permutation;
    Tensor<float,NDIM,BDIM> incident_matrix;
    Tensor<float,BDIM - 1,BDIM - 1> J;
    Tensor<float,SDIM,SDIM> Jx;
    Tensor<float,DDIM,DDIM> Jw;
    Tensor<float,PDIM,PDIM> Jy;
    Tensor<float,SDIM,DDIM> K;
    Tensor<float,SDIM,PDIM> Gx;
    Tensor<float,DDIM,PDIM> Gw;
    Tensor<float,SDIM,SDIM> Ix;
    Tensor<float,DDIM,DDIM> Iw;
    Tensor<float,SDIM,SDIM> Q;
    Tensor<float,SDIM,SDIM> D;
    Tensor<float,SDIM,SDIM> Ax;
    Tensor<float,SDIM,DDIM> Bx;
    Tensor<float,SDIM,PDIM> Cx;
    Tensor<float,DDIM,SDIM> Aw;
    Tensor<float,DDIM,DDIM> Bw;
    Tensor<float,DDIM,PDIM> Cw;
    Tensor<float,PDIM,SDIM> Ay;
    Tensor<float,PDIM,DDIM> By;
    Tensor<float,PDIM,PDIM> Cy;

    // Member Methods
    void setup();
    void computeMatrixJ();
    Tensor<float,DDIM> nonlinearFunction(const Tensor<float,DDIM>& w);
    Tensor<float,DDIM,DDIM> nonlinearJacobian(const Tensor<float,DDIM>& w);
    Tensor<float,DDIM> residualFunction(const Tensor<float,DDIM>& w);
    Tensor<float,DDIM,DDIM> residualJacobian(const Tensor<float,DDIM>& w);
    Tensor<float,DDIM> solveNonlinearity(const Tensor<float,DDIM>& w0);
};

// ====================================================================================================
// Class: PassiveFilter
// ====================================================================================================
class PassiveFilter final : public BaseModel
{
    public:
    // Constructor & Destructor Methods
    PassiveFilter();
    PassiveFilter(float fs,float tone = 0.5f);
    ~PassiveFilter();

    // Other Methods
    float processSample(float Vin) override;
    void setTone(float tone);
    float getTone();

    private:
    // Static Variables
    static constexpr int ORDER = 2;

    // Member Variables
    float tone;
    float R15;
    float R16;
    float R17;
    float VR2;
    float VR3;
    float C11;
    float C12;
    Tensor<float,ORDER + 1> x;
    Tensor<float,ORDER> y;
    Tensor<float,ORDER + 1> exponents;
    Tensor<float,ORDER + 1> k_B0;
    Tensor<float,ORDER + 1> k_B1;
    Tensor<float,ORDER + 1> k_B2;
    Tensor<float,ORDER + 1> k_A0;
    Tensor<float,ORDER + 1> k_A1;
    Tensor<float,ORDER + 1> k_A2;
    Tensor<float,ORDER + 1> b;
    Tensor<float,ORDER> a;

    // Member Methods
    void setup();
    void initializeCoeffs();
    void updateCoeffs();
};

// ====================================================================================================
// Class: BossDS1
// ====================================================================================================
class BossDS1 final : public BaseModel
{
    public:
    // Constructor & Destructor Methods
    BossDS1();
    BossDS1(float fs,float distortion = 0.5f,float tone = 0.5f,float level = 1.0f);
    ~BossDS1();

    // Other Methods
    float processSample(float Vin) override;
    void setDistortion(float distortion);
    void setTone(float tone);
    void setLevel(float level);
    float getDistortion();
    float getTone();
    float getLevel();

    private:
    // Member Variables
    float distortion;
    float tone;
    float level;
    TransistorAmplifier transistor_amplifier;
    OperationalAmplifier operational_amplifier;
    DiodeClipper diode_clipper;
    PassiveFilter passive_filter;
};

#endif