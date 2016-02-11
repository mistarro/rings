#ifndef RINGS__RMOD_HPP__
#define RINGS__RMOD_HPP__

#include <string>

#include "ZZ.hpp"

template <class R, unsigned ID>
class Rmod;

template <class R, unsigned ID>
std::string ToStringSpec(Rmod<R, ID> const &);

template <class R, unsigned ID>
Rmod<R, ID> operator-(Rmod<R, ID> const &);

template <class R, unsigned ID>
Rmod<R, ID> operator~(Rmod<R, ID> const &);

template <class R, unsigned ID>
class Rmod
{
public:
    Rmod() = default;
    Rmod(Rmod const &) = default;
    Rmod(Rmod &&) = default;

    template <typename S>
    Rmod(S const &);

    Rmod & operator=(Rmod const &) = default;
    Rmod & operator=(Rmod &&) = default;

    Rmod & operator+=(Rmod const &);
    Rmod & operator-=(Rmod const &);
    Rmod & operator*=(Rmod const &);

    friend Rmod<R, ID> operator-<R, ID>(Rmod<R, ID> const &);
    friend Rmod<R, ID> operator~<R, ID>(Rmod<R, ID> const &);

    static Rmod const & Zero();
    static Rmod const & One();

    bool operator==(Rmod const &) const;
    bool operator!=(Rmod const &) const;

    unsigned long int Hash() const;

    // output
    std::string ToString() const;
    friend std::string ToStringSpec<R, ID>(Rmod<R, ID> const & a);
    std::size_t StringLength() const;
    bool NeedParentheses() const;

    // performance
    static void Swap(Rmod &, Rmod &);

    // specific
    R const & Representative() const;

    static R const & Modulus();

    static void Init(R &&);

    // random
    static Rmod Random();

private:
    R value;

    static R modulus;

    static Rmod zero;
    static Rmod one;
};

template <class R, unsigned ID>
R Rmod<R, ID>::modulus(0);

template <class R, unsigned ID>
Rmod<R, ID> Rmod<R, ID>::zero(0);

template <class R, unsigned ID>
Rmod<R, ID> Rmod<R, ID>::one(1);

template <class R, unsigned ID>
template <typename S>
inline Rmod<R, ID>::Rmod(S const & a)
  : value(a)
{
    value %= modulus;
}

template <class R, unsigned ID>
inline Rmod<R, ID> & Rmod<R, ID>::operator+=(Rmod const & a)
{
    value += a.value;
    value %= modulus;
    return *this;
}

template <class R, unsigned ID>
inline Rmod<R, ID> & Rmod<R, ID>::operator-=(Rmod const & a)
{
    value -= a.value;
    value %= modulus;
    return *this;
}

template <class R, unsigned ID>
inline Rmod<R, ID> & Rmod<R, ID>::operator*=(Rmod const & a)
{
    value *= a.value;
    value %= modulus;
    return *this;
}

template <class R, unsigned ID>
inline Rmod<R, ID> const & Rmod<R, ID>::Zero()
{
    return zero;
}

template <class R, unsigned ID>
inline Rmod<R, ID> const & Rmod<R, ID>::One()
{
    return one;
}

template <class R, unsigned ID>
inline bool Rmod<R, ID>::operator==(Rmod const & a) const
{
    return (value == a.value);
}

template <class R, unsigned ID>
inline bool Rmod<R, ID>::operator!=(Rmod const & a) const
{
    return !(a == *this);
}

template <class R, unsigned ID>
inline Rmod<R, ID> operator-(Rmod<R, ID> const & a)
{
    return -a.value;
}

template <class R, unsigned ID>
inline Rmod<R, ID> operator~(Rmod<R, ID> const & a)
{
    return InvMod(a.value, a.modulus);
}

template <class R, unsigned ID>
inline std::string ToStringSpec(Rmod<R, ID> const & a)
{
    return a.value.ToString();
}

template <class R, unsigned ID>
inline std::string Rmod<R, ID>::ToString() const
{
	// every specialization may define its own ToStringSpec non-member function;
	// default is to return value.ToString() (see above);
    return ToStringSpec(*this);
}

template <class R, unsigned ID>
inline size_t Rmod<R, ID>::StringLength() const
{
    return value.StringLength();
}

template <class R, unsigned ID>
inline bool Rmod<R, ID>::NeedParentheses() const
{
    return value.NeedParentheses();
}

template <class R, unsigned ID>
inline unsigned long int Rmod<R, ID>::Hash() const
{
    return value.Hash() ^ (unsigned long int)ID;
}

template <class R, unsigned ID>
inline void Rmod<R, ID>::Swap(Rmod & a, Rmod & b)
{
    R::Swap(a.value, b.value);
}

template <class R, unsigned ID>
inline R const & Rmod<R, ID>::Representative() const
{
    return value;
}

template <class R, unsigned ID>
inline R const & Rmod<R, ID>::Modulus()
{
    return modulus;
}

template <class R, unsigned ID>
inline void Rmod<R, ID>::Init(R && m)
{
    modulus = std::forward<R>(m);
}

#endif // RINGS__RMOD_HPP__
