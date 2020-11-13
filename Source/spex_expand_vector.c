//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_expand_vector: double the space for a SPEX_vector object
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function expands a SPEX_vector by doubling its size. It will
 * initialize/allocate for the mpz entries if INIT_MPZ is true. Otherwise, this
 * aversion merely expands x and i and does not initialize/allocate the values!
 */

#include "spex_internal.h"

SPEX_info spex_expand_vector
(
    SPEX_vector* v, // the vector to be expanded
    bool INIT_MPZ   // indicate if the mpz entries should be initialized
)
{
    //--------------------------------------------------------------------------
    // double the size of v->x and v->i
    //--------------------------------------------------------------------------
    SPEX_info info;
    int64_t nzmax = v->nzmax ;

    bool okx, oki ;
    v->x = (mpz_t *)
        SPEX_realloc (2*nzmax, nzmax, sizeof (mpz_t), v->x, &okx) ;
    v->i = (int64_t *)
        SLIP_realloc (2*nzmax, nzmax, sizeof (int64_t), v->i, &oki) ;
    if (!oki || !okx)
    {
        return (SLIP_OUT_OF_MEMORY) ;
    }

    v->nzmax = 2*nzmax ;

    //--------------------------------------------------------------------------
    // set newly allocated mpz entries to NULL and initialize if required
    //--------------------------------------------------------------------------

    for (int64_t p = nzmax ; p < 2*nzmax ; p++)
    {
        SPEX_MPZ_SET_NULL (v->x[p]) ;
    }

    if (INIT_MPZ)
    {
        for (int64_t p = nzmax ; p < 2*nzmax ; p++)
        {
            SPEX_CHECK(SPEX_mpz_init (v->x[p])) ;
        }
    }

    return (SLIP_OK) ;
}
