#include "ZZ.hpp"

ZZ ZZ::zero(0);
ZZ ZZ::one(1);

unsigned long int ZZ::Hash() const
{
    mp_limb_t * limb = data[0]._mp_d;
    size_t size = mpz_size(data);
    mp_limb_t hash = size;
    for (int i = 0; i < size; ++i, ++limb)
        hash ^= *limb;
    return hash;
}

ZZ::RandState ZZ::state;
