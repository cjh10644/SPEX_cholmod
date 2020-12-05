//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_create_mpq_array: create a dense mpq array
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function creates an mpq array of size n, and initialize each of
 * the mpq_t entries.
 */

#include "spex_internal.h"

mpq_t* spex_create_mpq_array
(
    int64_t n            // size of the array
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (n <= 0) {return NULL;}

    //--------------------------------------------------------------------------

    // Malloc space
    mpq_t* x = (mpq_t*) SPEX_calloc(n, sizeof(mpq_t));
    if (!x) {return NULL;}
    for (int64_t i = 0; i < n; i++)
    {
        if (SPEX_mpq_init(x[i]) != SPEX_OK)
        {
            // Out of memory
            SPEX_MPQ_SET_NULL(x[i]);
            for (int64_t j = 0; j < n; j++)
            {
                if ( x[j] != NULL)
                {
                    SPEX_MPQ_CLEAR( x[j]);
                }
            }
            SPEX_FREE(x);
            return NULL;
        }
    }
    return x;
}
