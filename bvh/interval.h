/*
 Copyright (c) 2024 Andrew Astakhov <aa@contextmachine.ru>. All rights reserved.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE

*/


#ifndef INTERVAL_H
#define INTERVAL_H



#include <iostream>
#include <cmath>
#include <cstdlib>
#include <math.h>
namespace bvh{
#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

#ifndef MAYBE
#define MAYBE -1
#endif


#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#define FLOATING  0
#define INTERVAL  1

inline int INPUT_STYLE = INTERVAL;


typedef union du_sh {
    double d;
    unsigned short s[4];
} Double;
#ifdef WIN32
#define index_contain_exponent 3 // 0 for some UNIX platform 3 for some PC
#else
#define index_contain_exponent 0 // 0 for some UNIX platform 3 for some PC
#endif
#define MSW index_contain_exponent

inline du_sh RU_DB;
inline du_sh ROUND;

/*** Takashi's original version; correct but slow **********************/
/*** use the latest version of round() instead *************************/
inline double ulp(const double d)
{
    int expt;
    frexp(d, &expt);
    return ldexp(0.5, expt-52); // extract the ulp
}
/***********************************************************************/

class Interval {




public:
    double low;
    double upp;
    // Constructors
    Interval() : low(0.0), upp(0.0) {}
    Interval(const double value) : low(value), upp(value) {}
    Interval(const double low, const double upp) : low(low), upp(upp) {}
    Interval(const Interval& other) : low(other.low), upp(other.upp) {}

    // Destructor
    ~Interval() {}

    // Accessors
    double get_low() const { return low; }
    double get_upp() const { return upp; }
    void set_low(const double value) { low = value; }
    void set_upp(const double value) { upp = value; }

    // Methods
    double center() const {
        return (low + upp) / 2.0;
    }

    double range() const {
        return std::abs(upp - low);
    }

    double magnitude() const {
        return std::max(std::abs(low), std::abs(upp));
    }

    Interval lower() const {
        return Interval(low, low);
    }

    Interval upper() const {
        return Interval(upp, upp);
    }

    // Assignment operators
    Interval& operator=(const Interval& other) {
        if (this != &other) {
            low = other.low;
            upp = other.upp;
        }
        return *this;
    }

    Interval& operator=(const double value) {
        low = value;
        upp = value;
        return *this;
    }

    Interval& operator+=(const Interval& rhs) {
        low += rhs.low;
        upp += rhs.upp;
        return *this;
    }

    Interval& operator-=(const Interval& rhs) {
        low -= rhs.upp;
        upp -= rhs.low;
        return *this;
    }

    Interval& operator*=(const Interval& rhs) {
        double min_val = std::min({low * rhs.low, low * rhs.upp, upp * rhs.low, upp * rhs.upp});
        double max_val = std::max({low * rhs.low, low * rhs.upp, upp * rhs.low, upp * rhs.upp});
        low = min_val;
        upp = max_val;
        return *this;
    }

    Interval& operator/=(const Interval& rhs) {
        if (rhs.low <= 0.0 && rhs.upp >= 0.0) {
            throw std::runtime_error("Division by an interval containing zero is undefined.");
        }
        double min_val = std::min({low / rhs.low, low / rhs.upp, upp / rhs.low, upp / rhs.upp});
        double max_val = std::max({low / rhs.low, low / rhs.upp, upp / rhs.low, upp / rhs.upp});
        low = min_val;
        upp = max_val;
        return *this;
    }

    Interval& operator+=(const double value) {
        low += value;
        upp += value;
        return *this;
    }

    Interval& operator-=(const double value) {
        low -= value;
        upp -= value;
        return *this;
    }

    Interval& operator*=(const double value) {
        if (value >= 0) {
            low *= value;
            upp *= value;
        } else {
            double new_low = upp * value;
            double new_upp = low * value;
            low = new_low;
            upp = new_upp;
        }
        return *this;
    }

    Interval& operator/=(const double value) {
        if (value == 0.0) {
            throw std::runtime_error("Division by zero is undefined.");
        }
        if (value > 0) {
            low /= value;
            upp /= value;
        } else {
            double new_low = upp / value;
            double new_upp = low / value;
            low = new_low;
            upp = new_upp;
        }
        return *this;
    }

    // Arithmetic operators
    Interval operator+(const Interval& rhs) const {
        Interval result(*this);
        result += rhs;
        return result;
    }

    Interval operator-(const Interval& rhs) const {
        Interval result(*this);
        result -= rhs;
        return result;
    }

    Interval operator*(const Interval& rhs) const {
        Interval result(*this);
        result *= rhs;
        return result;
    }

    Interval operator/(const Interval& rhs) const {
        Interval result(*this);
        result /= rhs;
        return result;
    }

    Interval operator+(const double value) const {
        Interval result(*this);
        result += value;
        return result;
    }

    Interval operator-(const double value) const {
        Interval result(*this);
        result -= value;
        return result;
    }

    Interval operator*(const double value) const {
        Interval result(*this);
        result *= value;
        return result;
    }

    Interval operator/(const double value) const {
        Interval result(*this);
        result /= value;
        return result;
    }

    // int result enumeration
 

    // Conservative int methods
    inline int operator>(const Interval& rhs) const {
        if (low > rhs.upp)
            return TRUE;
        else if (upp <= rhs.low)
            return FALSE;
        else
            return MAYBE;
    }

    inline int operator<(const Interval& rhs) const {
        if (upp < rhs.low)
            return TRUE;
        else if (low >= rhs.upp)
            return FALSE;
        else
            return MAYBE;
    }

    inline int operator>=(const Interval& rhs) const {
        if (low >= rhs.upp)
            return TRUE;
        else if (upp < rhs.low)
            return FALSE;
        else
            return MAYBE;
    }

    inline int operator<=(const Interval& rhs) const {
        if (upp <= rhs.low)
            return TRUE;
        else if (low > rhs.upp)
            return FALSE;
        else
            return MAYBE;
    }

    // int with double
    inline int operator>(const double value) const {
        if (low > value)
            return TRUE;
        else if (upp <= value)
            return FALSE;
        else
            return MAYBE;
    }

    inline int operator<(const double value) const {
        if (upp < value)
            return TRUE;
        else if (low >= value)
            return FALSE;
        else
            return MAYBE;
    }

    inline int operator>=(const double value) const {
        if (low >= value)
            return TRUE;
        else if (upp < value)
            return FALSE;
        else
            return MAYBE;
    }

    inline int operator<=(const double value) const {
        if (upp <= value)
            return TRUE;
        else if (low > value)
            return FALSE;
        else
            return MAYBE;
    }

    // int with double on the left
    friend inline int operator>(const double value, const Interval& rhs) {
        if (value > rhs.upp)
            return TRUE;
        else if (value <= rhs.low)
            return FALSE;
        else
            return MAYBE;
    }

    friend inline int operator<(const double value, const Interval& rhs) {
        if (value < rhs.low)
            return TRUE;
        else if (value >= rhs.upp)
            return FALSE;
        else
            return MAYBE;
    }

    friend inline int operator>=(const double value, const Interval& rhs) {
        if (value >= rhs.upp)
            return TRUE;
        else if (value < rhs.low)
            return FALSE;
        else
            return MAYBE;
    }

    friend inline int operator<=(const double value, const Interval& rhs) {
        if (value <= rhs.low)
            return TRUE;
        else if (value > rhs.upp)
            return FALSE;
        else
            return MAYBE;
    }

    // Equality and inequality comparisons
    bool operator==(const Interval& rhs) const {
        return low == rhs.low && upp == rhs.upp;
    }

    bool operator!=(const Interval& rhs) const {
        return !(*this == rhs);
    }

    bool operator==(const double value) const {
        return low == value && upp == value;
    }

    bool operator!=(const double value) const {
        return !(*this == value);
    }

    // Other methods
    static Interval pow(const Interval& a, double n) {
        if (a.low < 0.0 && std::floor(n) != n) {
            throw std::runtime_error("Cannot raise a negative interval to a non-integer power.");
        }

        double min_val = std::pow(a.low, n);
        double max_val = std::pow(a.upp, n);

        if (n == int(n) && int(n) % 2 == 0 && a.low < 0.0 && a.upp > 0.0) {
            min_val = 0.0;
            max_val = std::max(std::pow(a.low, n), std::pow(a.upp, n));
        }

        if (min_val > max_val) std::swap(min_val, max_val);
        return Interval(min_val, max_val);
    }

    static Interval pow(const Interval& a, int n) {
        return Interval::pow(a, static_cast<double>(n));
    }

    static Interval sqrt(const Interval& a) {
        if (a.low < 0.0) {
            throw std::runtime_error("Square root of negative number is undefined.");
        }
        return Interval(std::sqrt(a.low), std::sqrt(a.upp));
    }

    static Interval merge(const Interval& a, const Interval& b) {
        double new_low = std::min(a.low, b.low);
        double new_upp = std::max(a.upp, b.upp);
        return Interval(new_low, new_upp);
    }
    /* this routine is to decide whether the num is in the interval. */
    bool in_interval(const double num)
    {

        if ( get_low() <= num && num <= get_upp())
            return TRUE;
        else
            return FALSE;
    }

    bool contain(const double num) const {
        if ( get_low() <= num && num <= get_upp())
            return TRUE;
        else
            return FALSE;

    }

    bool contain(const Interval& num) const {
        return num.low>=low && num.upp<=upp;


    }


    bool overlaps(const Interval& other) const {
        return !(low > other.upp || upp < other.low);
    }

    // Friend functions for stream operators
    friend std::istream& operator>>(std::istream& is, Interval& inv);
    friend std::ostream& operator<<(std::ostream& os, const Interval& inv);

    // Helper function to display int result
    static const char* comparison_to_string(int comp) {
        switch (comp) {
            case TRUE:
                return "TRUE";
            case FALSE:
                return "FALSE";
            case MAYBE:
                return "MAYBE";
            default:
                return "UNKNOWN";
        }
    }
};

// Non-member arithmetic operators with double on the left
inline Interval operator+(const double value, const Interval& rhs) {
    return rhs + value;
}

inline Interval operator-(const double value, const Interval& rhs) {
    Interval result(value - rhs.get_upp(), value - rhs.get_low());
    return result;
}

inline Interval operator*(const double value, const Interval& rhs) {
    return rhs * value;
}

inline Interval operator/(const double value, const Interval& rhs) {
    if (rhs.get_low() <= 0.0 && rhs.get_upp() >= 0.0) {
        throw std::runtime_error("Division by an interval containing zero is undefined.");
    }
    double min_val = std::min(value / rhs.get_low(), value / rhs.get_upp());
    double max_val = std::max(value / rhs.get_low(), value / rhs.get_upp());
    return Interval(min_val, max_val);
}

// Stream operators
inline std::ostream& operator<<(std::ostream& os, const Interval& inv) {
    os << "[" << inv.low << ", " << inv.upp << "]";
    return os;
}

inline std::istream& operator>>(std::istream& is, Interval& inv) {
    char ch;
    is >> ch; // Should be '['
    if (ch != '[') {
        is.setstate(std::ios::failbit);
        return is;
    }
    is >> inv.low;
    is >> ch; // Should be ','
    if (ch != ',') {
        is.setstate(std::ios::failbit);
        return is;
    }
    is >> inv.upp;
    is >> ch; // Should be ']'
    if (ch != ']') {
        is.setstate(std::ios::failbit);
        return is;
    }
    return is;
}

/*
double pow (double a, int n)
{
     double b = 1.0;
     if (n== 0 ) return b;


     for (int i=0; i< n ; i++)
     {
	  b = b * a;
     }
     return b;
}
*/
/*********************************************************************/
/*** return ulp of double precision number x (new version) ***********/
/*** correct and fast for normalized and denormalized ****************/

static unsigned short mask[16] = { /* bit masks for bits 0 - 15 */
  0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080,
  0x0100, 0x0200, 0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000};

inline double round(const double x)
{
  Double U,   /* ulp of x */
    X;        /* working copy of x */
  int bit,    /* position of bit e-1 in 16-bit word */
    e1,       /* biased exponent - 1 */
    word;     /* index of 16-bit word containing bit e-1 */

  X.d = x;
  X.s[MSW] &= 0x7ff0;              /* isolate exponent in 16-bit word */

  /* X.s[0] now holds the exponent in bits 14-4 */

  U.d = 0.0;                       /* initialize exponent and mantissa to 0 */

  if (X.s[MSW] > 0x0340)           /* ulp is normalized number */
    U.s[MSW] = X.s[MSW]-0x0340;    /*  set exponent to e-52 */

  /* the value 0x0340 is 52 left-shifted 4 bits, i.e. 0x0340 = 832 = 52<<4 */

  else {                           /* ulp is denormalized number */
    e1 = (X.s[MSW]>>4) - 1;        /*  biased exponent - 1 */
    word = e1>>4;                  /*  find 16-bit word containing bit e-1 */
    if (MSW == 0) word = 3 - word; /* compensate for word ordering */
    bit  = e1%16;                  /*  find the bit position in this word */

    U.s[word] |= mask[bit];        /*  set the bit to 1 */
  }

  return U.d;                      /* return ulp */
}


inline Interval merge(const Interval& a, const Interval& b)
  // Merge two Intervals
{
  double low;
  double upp;

  if (a.low < b.low)
    low = a.low;
  else
    low = b.low;

  if (a.upp > b.upp)
    upp = a.upp;
  else
    upp = b.upp;

  Interval inv; inv.low = low; inv.upp = upp; return(inv);
}

inline int if_overlap(const Interval& a, const Interval& b)
{
     if (a.low > b.upp) return FALSE;
     else if (a.upp < b.low) return FALSE;
     else return TRUE;

}


inline std::ostream& operator<<( std::ostream& os, Interval& inv)
{
   os << "[" << inv.low << ", " << inv.upp << "]";

//  os << "{" <<  ((inv.low + inv.upp) / 2.0) << ",}";
   return os;
}



  // --------------------
  // Arithmetic Operators
  // --------------------

  //---------------------
  // int Operators
  //---------------------


inline int comp_greater(Interval& a, Interval& b) // conservatively compare
  {
    if (a.get_low() > b.get_upp()) return TRUE;
    else if(a.get_upp() <= b.get_low()) return FALSE;
//    else return MAYBE;

    else if (a.get_upp() > b.get_upp()) return TRUE;
   else return FALSE;

  }


inline int comp_less(Interval& a, Interval& b) // conservatively compare
  {
    if (a.get_upp() < b.get_low()) return TRUE;
    else if(a.get_low() >= b.get_upp()) return FALSE;
//    else return MAYBE;
    else if (a.get_low() < b.get_low()) return TRUE;
    else return FALSE;
  }


    inline int operator > (double a, Interval& b)
  {
    if (a > b.upp) return TRUE;
    else if(a <= b.low) return FALSE;
    else return MAYBE;
  }

    inline int operator < (double a, Interval& b)
  {
       std::cout << " in < double, Interval " << std::endl;
    if (a < b.low) {std::cout << "out < d , I ,true" << std::endl;return TRUE;}
    else if(a >= b.upp) {std::cout << "out < d , I ,false" << std::endl;return FALSE;}
    else { std::cout << "out < d , I ,maybe" << std::endl;return MAYBE;}
  }

    inline int operator >= (double a, Interval& b)
  {
    if (a >= b.upp) return TRUE;
    else if(a <= b.low) return FALSE;
    else return MAYBE;
  }

    inline int operator <= (double a, Interval& b)
  {
    if (a <= b.low) return TRUE;
    else if(a >= b.upp) return FALSE;
    else return MAYBE;
  }


inline int operator == (double a, Interval& b)
  {
    return (a == b.low && a == b.upp);
  }

inline int operator != (double a, Interval& b)
  {
    return (a != b.low || a != b.upp);
  }
}

#endif //INTERVAL_H