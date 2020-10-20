//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_get_dense_Ak.c: build dense vector for column/row k of A
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to build dense mpz vector for column or 
// row k of A. Specifically, it is used to build Lk_dense_col and Uk_dense_row.


#include "spex_internal.h"

SPEX_info spex_get_dense_Ak
(
    mpz_t *Ak_dense_v,
    spex_vector *v,
)
{
    for (int64_t p = 0; p < v->nz; p++)
    {
        i = v->i[p];
        // TODO check if 0?
        SPEX_CHECK(SPEX_mpz_swap(Ak_dense_v[i], v->x[i]));
    }
}

