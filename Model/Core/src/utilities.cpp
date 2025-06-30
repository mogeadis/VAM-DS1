#include "utilities.h"

namespace utils
{
    // ====================================================================================================
    // Function: limit()
    // ====================================================================================================
    float limit(float value,float lower,float upper,float margin)
    {
        // Limit Value
        if(value <= lower)
        {
            value = lower + margin;
        }
        else if(value >= upper)
        {
            value = upper - margin;
        }

        // Return
        return value;
    }
}