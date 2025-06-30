#ifndef UTILITIES_H
#define UTILITIES_H

namespace utils
{
    // ====================================================================================================
    // Configuration Variables
    // ====================================================================================================
    
    // Newton-Raphson Parameters
    constexpr int ITERATIONS = 10;
    constexpr int SUBITERATIONS = 10;
    constexpr float THRESHOLD = 1e-6f;

    // Circuit Variables
    constexpr float POT_MARGIN = 1e-3f;
    constexpr float THERMAL_VOLTAGE = 25.85e-3f;

    // ====================================================================================================
    // Function: limit()
    // ====================================================================================================
    float limit(float value,float lower,float upper,float margin = 0);
}

#endif