#ifndef RINGS__COMMON_HPP__
#define RINGS__COMMON_HPP__

#include "ZZ.hpp"

/*
 *  Basic ring arithmetic using operators +=, -=, *=
 */

template <class R>
inline R operator+(R a, R const & b)
{
    return a += b;
}

template <class R>
inline R operator-(R a, R const & b)
{
    return a -= b;
}

template <class R>
inline R operator*(R a, R const & b)
{
    return a *= b;
}

// Euclidean rings: division with remainder
template <class R>
inline R operator/(R a, R const & b)
{
    return a /= b;
}

template <class R>
inline R operator%(R a, R const & b)
{
    return a %= b;
}

template <class R>
R operator^(R a, ZZ n)
{
    R res(R::One());

    if (n == 0)
        return res;

    if (n < 0)
    {
        n = -n;
        a = ~a;
    }

    // compute a^n, now n > 0
    for (SmallNat i = n.BitLength(); i-- > 0;)
    {
        res *= res;
        if (n[i])
            res *= a;
    }

    return res;
}

#endif // RINGS__COMMON_HPP__
