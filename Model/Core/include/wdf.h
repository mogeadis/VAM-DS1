#ifndef WDF_H
#define WDF_H

#include <array>

#include <Fastor/Fastor.h>

using namespace Fastor;

namespace wdf
{
    // ====================================================================================================
    // Class: BasePort
    // ====================================================================================================
    class BasePort
    {
        public:
        // Constructor & Destructor Methods
        BasePort();
        BasePort(float Rp);
        virtual ~BasePort();

        // Other Methods
        virtual float getReflectedWave() = 0;
        virtual void setIncidentWave(float a) final;
        virtual void setPortResistance(float Rp) final;
        virtual float getPortResistance() final;
        virtual float waveToVoltage() final;
        virtual float waveToCurrent() final;

        protected:
        // Member Variables
        float Rp;
        float a;
        float b;

        private:
        // Member Methods
        void setup();
    };

    // ====================================================================================================
    // Class: Resistor
    // ====================================================================================================
    class Resistor final : public BasePort
    {
        public:
        // Constructor & Destructor Methods
        Resistor();
        Resistor(float R);
        ~Resistor();

        // Other Methods
        float getReflectedWave() override;
    };

    // ====================================================================================================
    // Class: Capacitor
    // ====================================================================================================
    class Capacitor final : public BasePort
    {
        public:
        // Constructor & Destructor Methods
        Capacitor();
        Capacitor(float C,float fs = 44100.0f);
        ~Capacitor();

        // Other Methods
        float getReflectedWave() override;
    };

    // ====================================================================================================
    // Class: ResistiveVoltageSource
    // ====================================================================================================
    class ResistiveVoltageSource final : public BasePort
    {
        public:
        // Constructor & Destructor Methods
        ResistiveVoltageSource();
        ResistiveVoltageSource(float Rs);
        ~ResistiveVoltageSource();

        // Other Methods
        float getReflectedWave() override;
        void setVoltage(float Vs);

        private:
        // Member Variables
        float Vs;

        // Member Methods
        void setup();
    };

    // ====================================================================================================
    // Class: ResistiveCurrentSource
    // ====================================================================================================
    class ResistiveCurrentSource final : public BasePort
    {
        public:
        // Constructor & Destructor Methods
        ResistiveCurrentSource();
        ResistiveCurrentSource(float Rs);
        ~ResistiveCurrentSource();

        // Other Methods
        float getReflectedWave() override;
        void setCurrent(float Is);

        private:
        // Member Variables
        float Is;

        // Member Methods
        void setup();
    };

    // ====================================================================================================
    // Class: BaseAdaptorRoot
    // ====================================================================================================

    template <typename...Ports>
    class BaseAdaptorRoot
    {
        public:
        // Constructor & Destructor Methods
        BaseAdaptorRoot();
        BaseAdaptorRoot(Ports&...ports);
        virtual ~BaseAdaptorRoot();

        // Other Methods
        virtual void computeScatteringMatrix() = 0;
        virtual void collectIncidentWaves() final;
        virtual void propagateReflectedWaves() final;

        protected:
        // Member Variables
        std::vector<BasePort*> ports; 
        int N;
        Fastor::Tensor<float,sizeof...(Ports)> a;
        Fastor::Tensor<float,sizeof...(Ports)> b;
        Fastor::Tensor<float,sizeof...(Ports),sizeof...(Ports)> scattering_matrix;
    };

    // ====================================================================================================
    // Class: SeriesAdaptorRoot
    // ====================================================================================================

    template <typename... Ports>
    class SeriesAdaptorRoot : public BaseAdaptorRoot<Ports...>
    {
        public:
        // Constructor & Destructor Methods
        SeriesAdaptorRoot();
        SeriesAdaptorRoot(Ports&...ports);
        ~SeriesAdaptorRoot();

        // Other Methods
        void computeScatteringMatrix() override;
    };

    // ====================================================================================================
    // Class: ParallelAdaptorRoot
    // ====================================================================================================

    template <typename... Ports>
    class ParallelAdaptorRoot : public BaseAdaptorRoot<Ports...>
    {
        public:
        // Constructor & Destructor Methods
        ParallelAdaptorRoot();
        ParallelAdaptorRoot(Ports&...ports);
        ~ParallelAdaptorRoot();

        // Other Methods
        void computeScatteringMatrix() override;
    };
}

#include "wdf.tpp"

#endif