//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_expand_vector: expand the space for a SPEX_vector object to
// given new size.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function expands a SPEX_vector to given size. It will
 * initialize/allocate for the mpz entries.
 */

#include "spex_internal.h"

SPEX_info spex_expand_vector
(
    SPEX_vector* v, // the vector to be expanded
    const int64_t new_size// desired new size for v
)
{
    //--------------------------------------------------------------------------
    // expand the size of v->x and v->i to new_size
    //--------------------------------------------------------------------------
    SPEX_info info;
    int64_t old_size = v->nzmax ;

    bool okx, oki ;
    v->x = (mpz_t *)
        SPEX_realloc (new_size, old_size, sizeof (mpz_t), v->x, &okx) ;
    v->i = (int64_t *)
        SPEX_realloc (new_size, old_size, sizeof (int64_t), v->i, &oki) ;
    if (!oki || !okx)
    {
        return (SPEX_OUT_OF_MEMORY) ;
    }

    v->nzmax = new_size ;

    //--------------------------------------------------------------------------
    // set newly allocated mpz entries to NULL and initialize if required
    //--------------------------------------------------------------------------

    for (int64_t p = old_size ; p < new_size ; p++)
    {
        SPEX_MPZ_SET_NULL (v->x[p]) ;
    }

    for (int64_t p = old_size ; p < new_size ; p++)
    {
        SPEX_CHECK(SPEX_mpz_init (v->x[p])) ;
    }

    return (SPEX_OK) ;
}
