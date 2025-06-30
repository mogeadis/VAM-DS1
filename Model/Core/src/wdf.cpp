#include "wdf.h"

namespace wdf
{
    // ====================================================================================================
    // Class: BasePort
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    BasePort::BasePort()
    {
        this->Rp = 1.0f;
        this->setup();
    }
    BasePort::BasePort(float Rp)
    {
        this->Rp = Rp;
        this->setup();
    }
    BasePort::~BasePort() = default;

    // [Public] Other Methods
    void BasePort::setIncidentWave(float a)
    {
        this->a = a;
    }
    void BasePort::setPortResistance(float Rp)
    {
        this->Rp = Rp;
    }
    float BasePort::getPortResistance()
    {
        return this->Rp;
    }
    float BasePort::waveToVoltage()
    {
        // Calculate Voltage
        float voltage = (this->a + this->b)/2;

        // Return
        return voltage;
    }
    float BasePort::waveToCurrent()
    {
        // Calculate Current
        float current = (this->a - this->b)/(2*this->Rp);
        
        // Return
        return current;
    }

    // [Private] Member Methods
    void BasePort::setup()
    {
        this->a = 0.0f;
        this->b = 0.0f;
    }

    // ====================================================================================================
    // Class: Resistor
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    Resistor::Resistor() : BasePort() {}
    Resistor::Resistor(float R) : BasePort(R) {}
    Resistor::~Resistor() = default;

    // [Public] Other Methods
    float Resistor::getReflectedWave()
    {
        this->b = 0;
        return this->b;
    }

    // ====================================================================================================
    // Class: Capacitor
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    Capacitor::Capacitor() : BasePort() {}
    Capacitor::Capacitor(float C,float fs) : BasePort(1/(2*fs*C)) {}
    Capacitor::~Capacitor() = default;

    // [Public] Other Methods
    float Capacitor::getReflectedWave()
    {
        this->b = this->a;
        return this->b;
    }

    // ====================================================================================================
    // Class: ResistiveVoltageSource
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    ResistiveVoltageSource::ResistiveVoltageSource() : BasePort()
    {
        this->setup();
    }
    ResistiveVoltageSource::ResistiveVoltageSource(float Rs) : BasePort(Rs)
    {
        this->setup();
    }
    ResistiveVoltageSource::~ResistiveVoltageSource() = default;

    // [Public] Other Methods
    float ResistiveVoltageSource::getReflectedWave()
    {
        this->b = this->Vs;
        return this->b;
    }
    void ResistiveVoltageSource::setVoltage(float Vs)
    {
        this->Vs = Vs;
    }

    // [Private] Member Methods
    void ResistiveVoltageSource::setup()
    {
        this->Vs = 0;
    }

    // ====================================================================================================
    // Class: ResistiveCurrentSource
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    ResistiveCurrentSource::ResistiveCurrentSource() : BasePort()
    {
        this->setup();
    }
    ResistiveCurrentSource::ResistiveCurrentSource(float Rs) : BasePort(Rs)
    {
        this->setup();
    }
    ResistiveCurrentSource::~ResistiveCurrentSource() = default;

    // [Public] Other Methods
    float ResistiveCurrentSource::getReflectedWave()
    {
        this->b = this->Is*this->Rp;
        return this->b;
    }
    void ResistiveCurrentSource::setCurrent(float Is)
    {
        this->Is = Is;
    }
    
    // [Private] Member Methods
    void ResistiveCurrentSource::setup()
    {
        this->Is = 0;
    }
}