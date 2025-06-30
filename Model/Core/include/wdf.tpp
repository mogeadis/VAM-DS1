#ifndef WDF_TPP
#define WDF_TPP

#include "wdf.h"

namespace wdf
{
    // ====================================================================================================
    // Class: BaseAdaptorRoot
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    template <typename...Ports>
    BaseAdaptorRoot<Ports...>::BaseAdaptorRoot()
    {
        this->N = 0;
    }
    template <typename...Ports>
    BaseAdaptorRoot<Ports...>::BaseAdaptorRoot(Ports&...ports)
    {
        (this->ports.push_back(&ports),...);
        this->N = static_cast<int>(this->ports.size());
        this->a.zeros();
        this->b.zeros();
        this->scattering_matrix.zeros();
    }
    template <typename...Ports>
    BaseAdaptorRoot<Ports...>::~BaseAdaptorRoot() = default;

    // [Public] Other Methods
    template <typename...Ports>
    void BaseAdaptorRoot<Ports...>::collectIncidentWaves()
    {
        for(int index = 0; index < this->N; index++)
        {
            this->a(index) = this->ports[index]->getReflectedWave();
        }
    }
    template <typename...Ports>
    void BaseAdaptorRoot<Ports...>::propagateReflectedWaves()
    {
        this->b = this->scattering_matrix % this->a;
        for(int index = 0; index < this->N; index++)
        {
            this->ports[index]->setIncidentWave(this->b(index));
        }
    }

    // ====================================================================================================
    // Class: SeriesAdaptorRoot
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    template <typename...Ports>
    SeriesAdaptorRoot<Ports...>::SeriesAdaptorRoot() : BaseAdaptorRoot<Ports...>()
    {
        this->computeScatteringMatrix();
    }
    template <typename...Ports>
    SeriesAdaptorRoot<Ports...>::SeriesAdaptorRoot(Ports&...ports) : BaseAdaptorRoot<Ports...>(ports...)
    {
        this->computeScatteringMatrix();
    }
    template <typename...Ports>
    SeriesAdaptorRoot<Ports...>::~SeriesAdaptorRoot() = default;

    // [Public] Other Methods
    template <typename...Ports>
    void SeriesAdaptorRoot<Ports...>::computeScatteringMatrix()
    {
        // Precompute Sum
        float sum = 0;
        for(auto port : this->ports)
        {
            sum += port->getPortResistance();
        }

        // Fill Scattering Matrix
        for(int row = 0; row < this->N; row++)
        {
            float coeff = -2*this->ports[row]->getPortResistance()/sum;
            for(int col = 0; col < this->N; col++)
            {
                this->scattering_matrix(row,col) = (row == col) ? coeff + 1 : coeff;
            }
        }
    }

    // ====================================================================================================
    // Class: ParallelAdaptorRoot
    // ====================================================================================================

    // [Public] Constructor & Destructor Methods
    template <typename...Ports>
    ParallelAdaptorRoot<Ports...>::ParallelAdaptorRoot() : BaseAdaptorRoot<Ports...>()
    {
        this->computeScatteringMatrix();
    }
    template <typename...Ports>
    ParallelAdaptorRoot<Ports...>::ParallelAdaptorRoot(Ports&...ports) : BaseAdaptorRoot<Ports...>(ports...)
    {
        this->computeScatteringMatrix();
    }
    template <typename...Ports>
    ParallelAdaptorRoot<Ports...>::~ParallelAdaptorRoot() = default;

    // [Public] Other Methods
    template <typename...Ports>
    void ParallelAdaptorRoot<Ports...>::computeScatteringMatrix()
    {
        // Precompute Sum
        float sum = 0;
        for(auto port : this->ports)
        {
            sum += 1/port->getPortResistance();
        }

        // Fill Scattering Matrix
        for(int col = 0; col < this->N; col++)
        {
            float coeff = 2*(1/this->ports[col]->getPortResistance())/sum;
            for(int row = 0; row < this->N; row++)
            {
                this->scattering_matrix(row,col) = (row == col) ? coeff - 1 : coeff;
            }
        }
    }
}

#endif