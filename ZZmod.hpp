#ifndef RINGS__ZZMOD_HPP__
#define RINGS__ZZMOD_HPP__

#include "ZZ.hpp"
#include "Rmod.hpp"
#include "common.hpp"

template <unsigned ID>
using ZZmod = Rmod<ZZ, ID>;

template <unsigned ID>
inline ZZmod<ID> operator-(ZZmod<ID> const & a)
{
    return ZZmod<ID>::Modulus() - a.Representative();
}

template <unsigned ID>
inline ZZmod<ID> operator^(ZZmod<ID> const & a, ZZ const & e)
{
    return PowMod(a.Representative(), e, ZZmod<ID>::Modulus());
}

#endif // RINGS__ZZMOD_HPP__
