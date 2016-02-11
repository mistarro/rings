#ifndef RINGS__ZZ_HPP__
#define RINGS__ZZ_HPP__

#include <string>
#include <ctime>
#include <gmp.h>

long int typedef SmallInt;
unsigned long int typedef SmallNat;

class ZZ
{
public:
    ZZ() { mpz_init(data); }
    ZZ(ZZ const & a) { mpz_init_set(data, a.data); }
    ZZ(ZZ && a) { mpz_init(data); mpz_swap(data, a.data); }
    ZZ(int a) { mpz_init_set_si(data, a); }
    explicit ZZ(SmallInt a) { mpz_init_set_si(data, a); }
    explicit ZZ(SmallNat a) { mpz_init_set_ui(data, a); }

    ~ZZ() { mpz_clear(data); }

    ZZ & operator=(ZZ const & a) { mpz_set(data, a.data); return *this; } // self-assignment ok
    ZZ & operator=(ZZ && a) { mpz_swap(data, a.data); return *this; }

    ZZ & operator+=(ZZ const & a) { mpz_add(data, data, a.data); return *this; }
    ZZ & operator-=(ZZ const & a) { mpz_sub(data, data, a.data); return *this; }
    ZZ & operator*=(ZZ const & a) { mpz_mul(data, data, a.data); return *this; }

    friend ZZ operator-(ZZ const & a) { ZZ res; mpz_neg(res.data, a.data); return res; }
    friend ZZ operator~(ZZ const & a) { return (a == 1 || a == -1) ? a : zero; }

    static ZZ const & Zero() { return zero; }
    static ZZ const & One() { return one; }

    bool operator==(ZZ const & a) const { return mpz_cmp(data, a.data) == 0; }
    bool operator!=(ZZ const & a) const { return mpz_cmp(data, a.data) != 0; }

    // output
    std::string ToString(int base = 10) const
    {
        char cstr[mpz_sizeinbase(data, base) + 1];
        mpz_get_str(cstr, base, data);
        std::string str(cstr);
        return str;
    }

    std::size_t StringLength(int base = 10) const { return mpz_sizeinbase(data, base); }
    bool NeedParentheses() const { return false; }

    unsigned long int Hash() const;

    // performance
    static void Swap(ZZ & a, ZZ & b) { mpz_swap(a.data, b.data); }

    // ZZ specific
    ZZ & operator/=(ZZ const & a) { mpz_fdiv_q(data, data, a.data); return *this; }
    ZZ & operator%=(ZZ const & a) { mpz_fdiv_r(data, data, a.data); return *this; }

    ZZ operator^(unsigned long int e) const { ZZ res; mpz_pow_ui(res.data, data, e); return res; }

    ZZ & operator++() { mpz_add_ui(data, data, 1); return *this; }
    ZZ & operator--() { mpz_sub_ui(data, data, 1); return *this; }

    bool operator>(ZZ const & a) const { return mpz_cmp(data, a.data) > 0; }
    bool operator<(ZZ const & a) const { return mpz_cmp(data, a.data) < 0; }
    bool operator>=(ZZ const & a) const { return mpz_cmp(data, a.data) >= 0; }
    bool operator<=(ZZ const & a) const { return mpz_cmp(data, a.data) <= 0; }

    bool operator[](SmallNat a) const { return mpz_tstbit(data, (mp_bitcnt_t)a); }

    SmallNat BitLength() const { return (SmallNat)mpz_sizeinbase(data, 2); }

    ZZ(const char * str, int base = 0) { mpz_init(data); mpz_set_str(data, str, base); }

    // misc conversions
    SmallNat Reduce(SmallNat a) const { return mpz_fdiv_ui(data, a);  }
    SmallInt ToSmallInt() const { return mpz_get_si(data); }
    SmallNat ToSmallNat() const { return mpz_get_ui(data); }

    // random
    static ZZ Random(SmallNat bits) { ZZ res; mpz_urandomb(res.data, state.s, bits); return res; }
    static ZZ Random(ZZ const & n) { ZZ res; mpz_urandomm(res.data, state.s, n.data); return res; }

    // friends
    friend void DivMod(ZZ const &, ZZ const &, ZZ &, ZZ &);
    friend ZZ Gcd(ZZ const &, ZZ const &);
    friend ZZ ExtGcd(ZZ const &, ZZ const &, ZZ &, ZZ &);
    friend ZZ InvMod(ZZ const &, ZZ const &);
    friend ZZ PowMod(ZZ const &, ZZ const &, ZZ const &);
    friend bool Root(ZZ const &, SmallNat, ZZ &);
    friend bool IsPrimeMR(ZZ const &, int);
    friend bool IsPerfectPower(ZZ const &);

    friend SmallNat operator/(SmallNat, ZZ const &);
    friend SmallInt operator/(SmallInt, ZZ const &);

private:
    mpz_t data;

    static ZZ zero;
    static ZZ one;

    // random generation
    class RandState
    {
    public:
        RandState() { gmp_randinit_default(s); gmp_randseed_ui(s, static_cast<SmallNat>(std::time(0))); }
        ~RandState() { gmp_randclear(s); }
        friend class ZZ;
    private:
        gmp_randstate_t s;
    };

    static RandState state;
};

// Auxiliary arithmetic functions

inline void DivMod(ZZ const & a, ZZ const & b, ZZ & q, ZZ & r)
{
    mpz_fdiv_qr(q.data, r.data, a.data, b.data);
}

inline ZZ Gcd(ZZ const & a, ZZ const & b)
{
    ZZ res;
    mpz_gcd(res.data, a.data, b.data);
    return res;
}

inline ZZ ExtGcd(ZZ const & a, ZZ const & b, ZZ & s, ZZ & t)
{
    ZZ res;
    mpz_gcdext(res.data, s.data, t.data, a.data, b.data);
    return res;
}

inline ZZ InvMod(ZZ const & a, ZZ const & n)
{
    ZZ res;
    if (!mpz_invert(res.data, a.data, n.data))
        return ZZ::zero;
    return res;
}

inline ZZ PowMod(ZZ const & a, ZZ const & e, ZZ const & n)
{
    if (e < 0 && Gcd(a, n) != 1)
        return ZZ::zero;
    ZZ res;
    mpz_powm(res.data, a.data, e.data, n.data);
    return res;
}

// Number Theoretic functions

// Miller-Rabin randomized primality test
inline bool IsPrimeMR(ZZ const & n, int reps = 7)
{
    return (mpz_probab_prime_p(n.data, reps) > 0);
}

inline bool IsPerfectPower(ZZ const & n)
{
    return (n != 0 && n != 1 && mpz_perfect_power_p(n.data));
}

// Other functions

inline bool Root(ZZ const & n, SmallNat k, ZZ & r)
{
    return mpz_root(r.data, n.data, k);
}

inline SmallNat operator/(SmallNat a, ZZ const & b)
{
    mpz_t aux;
    mpz_init_set_ui(aux, a);
    mpz_fdiv_q(aux, aux, b.data);
    // result fits into SmallNat
    return mpz_get_ui(aux);
}

inline SmallInt operator/(SmallInt a, ZZ const & b)
{
    mpz_t aux;
    mpz_init_set_si(aux, a);
    mpz_fdiv_q(aux, aux, b.data);
    // result fits into SmallInt
    return mpz_get_si(aux);
}

#include "common.hpp"

#endif // RINGS__ZZ_HPP__
