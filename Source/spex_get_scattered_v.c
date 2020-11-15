//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_get_scattered_v.c: build scattered vector for given sparse
// vector
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to build scattered mpz vector for column or
// row k of A, or the inserted column. This function eliminates explicit
// 0 if specified. If keep_v is false, this function swap the mpz values, and
// thus the original vector will become all zeros. Otherwise, mpz_set will be
// used to keep the original mpz values.


#include "spex_internal.h"

SPEX_info spex_get_scattered_v
(
    spex_scattered_vector *sv,   // output vector in scattered form
    const spex_vector *v,        // the vector in compressed form, whose
                                 // max index is n
    const int64_t n,             // number of entries in v
    const bool eliminate_zero,   // indicate if explicit zero should be elimated
    const bool keep_v            // indicate if the mpz values should be kept
)
{
    if (!v)
    {
        return SPEX_INCORRECT_INPUT;
    }

    sv = spex_create_scattered_vector(n);
    if (!sv)
    {
        return SPEX_OUT_OF_MEMORY;
    }

    SPEX_info info;
    int64_t p, i;
    /* TODO combine spex_find_next_nz with this function?
    int sgn;
    p = 0;
    while (p < v->nz)
    {
        i = v->i[p];
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, v->x[i]));
        if (sgn == 0 && eliminate_zero)
        {
            v->nz--;
            v->i[p] = v->i[v->nz];
        }
        else
        {
            if (sgn != 0 && k != perm_inv[i] && perm_inv[i] < inext)
            {
                inext = i;
            }
            if (keep_v)
            {
                SPEX_CHECK(SPEX_mpz_set(sv->x[i], v->x[i]));
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_swap(sv->x[i], v->x[i]));
            }
            sv->i[p] = v->i[p];
            p++;
        }
    }
    sv->nz = v->nz;
    */

    if (eliminate_zero)
    {
        int sgn;
        p = 0;
        while (p < v->nz)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, v->x[p]));
            if (sgn != 0)
            {
                i = v->i[p];
                if (keep_v)
                {
                    SPEX_CHECK(SPEX_mpz_set(sv->x[i], v->x[p]));
                }
                else
                {
                    SPEX_CHECK(SPEX_mpz_swap(sv->x[i], v->x[p]));
                }
                sv->i[p] = v->i[p];
                p++;
            }
            else
            {
                v->nz--;
                SPEX_CHECK(SPEX_mpz_swap(v->x[p], v->x[v->nz]));
                v->i[p] = v->i[v->nz];
            }
        }
        sv->nz = v->nz;
    }
    else
    {
        for (p = 0; p < v->nz; p++)
        {
            i = v->i[p];
            if (keep_v)
            {
                SPEX_CHECK(SPEX_mpz_set(sv->x[i], v->x[p]));
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_swap(sv->x[i], v->x[p]));
            }
            sv->i[p] = v->i[p];
        }
        sv->nz = v->nz;
    }
    return SPEX_OK;
}

