#ifndef VEC_H
#define VEC_H
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <vector>
#include <cstdlib>    // for std::abort
#include <type_traits>
#include <limits>
namespace bvh {

/******************************************************************************************
 *  detail::almost_equal:  Epsilon-based floating comparison or direct equality.
 *****************************************************************************************/
namespace detail {
template <typename T>
constexpr bool almost_equal(T a, T b, T eps = std::numeric_limits<T>::epsilon()) {
    if constexpr (std::is_floating_point_v<T>) {
        return (std::fabs(a - b) <= eps);
    } else {
        return (a == b);
    }
}
} // end namespace detail

/******************************************************************************************
 *  Primary Template (N-D) - Minimal fallback for arbitrary dimensions
 *****************************************************************************************/
template <typename T, std::size_t N>

class vec {
public:
    std::array<T, N> data_{};
    using value_type = T;

    using size_type = std::size_t;

    static constexpr std::size_t dim = N;
    constexpr vec()  {
        for(std::size_t i = 0; i < N; ++i) {
            data_[i] = T{};
        }
    }
    constexpr explicit vec(T val)  {
        for(std::size_t i = 0; i < N; ++i) {
            data_[i] = val;
        }
    }
    constexpr explicit vec(const std::array<T, N>& arr)  : data_{arr} {}

    // Construct from N arguments
    template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) == N)>>
    constexpr vec(Args... args)  : data_{static_cast<T>(args)...} {}

    constexpr T& operator[](std::size_t i) {
        if (i >= N) {
            throw std::out_of_range("vec index out of range");
        }
        return data_[i];
    }
    constexpr const T& operator[](std::size_t i) const {
        if (i >= N) {
            throw std::out_of_range("vec index out of range");
        }
        return data_[i];
    }

    constexpr bool operator==(const vec<T, N>& other) const  {
        for(std::size_t i = 0; i < N; ++i) {
            if(!detail::almost_equal(data_[i], other.data_[i])) {
                return false;
            }
        }
        return true;
    }
    constexpr bool operator!=(const vec<T, N>& other) const  {
        return !(*this == other);
    }

    constexpr vec<T, N> operator-() const  {
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            result.data_[i] = -data_[i];
        }
        return result;
    }
    constexpr vec<T, N> operator+(const vec<T, N>& rhs) const  {
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            result.data_[i] = data_[i] + rhs.data_[i];
        }
        return result;
    }
    constexpr vec<T, N>& operator+=(const vec<T, N>& rhs)  {
        for(std::size_t i = 0; i < N; ++i) {
            data_[i] += rhs.data_[i];
        }
        return *this;
    }
    constexpr vec<T, N> operator-(const vec<T, N>& rhs) const  {
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            result.data_[i] = data_[i] - rhs.data_[i];
        }
        return result;
    }
    constexpr vec<T, N>& operator-=(const vec<T, N>& rhs)  {
        for(std::size_t i = 0; i < N; ++i) {
            data_[i] -= rhs.data_[i];
        }
        return *this;
    }
    constexpr vec<T, N> operator*(T scalar) const  {
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            result.data_[i] = data_[i] * scalar;
        }
        return result;
    }
    constexpr vec<T, N>& operator*=(T scalar)  {
        for(std::size_t i = 0; i < N; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }
    constexpr vec<T, N> operator/(T val) const {
        if (detail::almost_equal(val, T{})) {
            throw std::invalid_argument("Division by zero");
        }
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            result.data_[i] = data_[i] / val;
        }
        return result;
    }
    constexpr vec<T, N>& operator/=(T val) {
        if (detail::almost_equal(val, T{})) {
            throw std::invalid_argument("Division by zero");
        }
        for(std::size_t i = 0; i < N; ++i) {
            data_[i] /= val;
        }
        return *this;
    }
    // Component-wise multiply
    constexpr vec<T, N> operator*(const vec<T, N>& rhs) const  {
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            result.data_[i] = data_[i] * rhs.data_[i];
        }
        return result;
    }
    // Component-wise divide
    constexpr vec<T, N> operator/(const vec<T, N>& rhs) const {
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            if (detail::almost_equal(rhs.data_[i], T{})) {
                throw std::invalid_argument("Division by zero in component-wise division");
            }
            result.data_[i] = data_[i] / rhs.data_[i];
        }
        return result;
    }

    // Utilities
    constexpr void set(T val)  {
        for(std::size_t i = 0; i < N; ++i) {
            data_[i] = val;
        }
    }
    template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) == N)>>
    constexpr void set(Args... args)  {
        data_ = {static_cast<T>(args)...};
    }
    constexpr void set(const std::array<T, N>& arr)  {
        data_ = arr;
    }
    constexpr T sqLength() const  {
        T sum = T{};
        for(std::size_t i = 0; i < N; ++i) {
            sum += data_[i]*data_[i];
        }
        return sum;
    }
    constexpr T length() const {
        return static_cast<T>(std::sqrt(sqLength()));
    }
    constexpr void unitize() {
        T l = length();
        if(detail::almost_equal(l, T{})) {
            throw std::invalid_argument("Cannot normalize zero-length vector");
        }
        for(std::size_t i = 0; i < N; ++i) {
            data_[i] /= l;
        }
    }

    /// Non-mutating normalization
    constexpr vec<T, N> unit() const {
        T l = length();
        if(detail::almost_equal(l, T{})) {
            throw std::invalid_argument("Cannot normalize zero-length vector");
        }
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            result.data_[i] = data_[i] / l;
        }
        return result;
    }

    /// Dot product with another vector
    constexpr T dot(const vec<T, N>& other) const  {
        T sum = T{};
        for(std::size_t i = 0; i < N; ++i) {
            sum += data_[i]*other.data_[i];
        }
        return sum;
    }

    /// Distance to another vector
    constexpr T distance(const vec<T, N>& other) const {
        return (other - *this).length();
    }

    /// Distance squared
    constexpr T distanceSq(const vec<T, N>& other) const  {
        return (other - *this).sqLength();
    }

    /// Projection length of this onto 'other'
    /// result = (this·other) / (other·other)
    constexpr T projection(const vec<T, N>& other) const {
        T dd = other.dot(other);
        if(detail::almost_equal(dd, T{})) {
            throw std::invalid_argument("Cannot project onto zero-length vector");
        }
        return dot(other)/dd;
    }

    /// vec projection of this onto 'other':
    /// result = other * [ (this·other) / (other·other) ]
    constexpr vec<T, N> project(const vec<T, N>& other) const {
        T dd = other.dot(other);
        if(detail::almost_equal(dd, T{})) {
            throw std::invalid_argument("Cannot project onto zero-length vector");
        }
        return other*(dot(other)/dd);
    }
    constexpr bool collinear(const vec<T, N>& other) const {
        T l1 = length();
        T l2 = other.length();
        if(detail::almost_equal(l1, T{}) || detail::almost_equal(l2, T{})) {
            return false;
        }
        T dotVal = T{};
        for(std::size_t i=0; i<N; ++i) {
            dotVal += (data_[i]/l1)*(other.data_[i]/l2);
        }
        return detail::almost_equal(std::fabs(dotVal), T(1.0));
    }

    /*********************
     *  Streaming
     *********************/
    friend std::ostream& operator<<(std::ostream& os, const vec<T, N>& vec) {
        os << "[";
        for(std::size_t i = 0; i < N; ++i) {
            os << vec.data_[i];
            if(i+1 < N) os << ",";
        }
        os << "]";
        return os;
    }
    vec abs() {
        vec<T, N> result;
        for(std::size_t i = 0; i < N; ++i) {
            result[i] = std::abs(data_[i]);
        }
        return result;
    }
    T max(const T val ) const{
        T max_val=val;
        for (int i = 0; i < N; ++i) {
            max_val=std::max(max_val,data_[i]);

        }
        return max_val;
    }
    vec<T,N> max(const vec<T,N> &val) const {

        vec<T,N> max_val;

        for (int i = 0; i < N; ++i) {
            max_val[i]=std::max(val[i],data_[i]);

        }
        return max_val;

    }
    T min(const T val ) const{
        T min_val=val;
        for (int i = 0; i < N; ++i) {
            min_val=std::min(min_val,data_[i]);

        }
        return min_val;

    }
    vec<T,N> min(const vec<T,N> &val) const {

        vec<T,N> min_val;

        for (int i = 0; i < N; ++i) {
            min_val[i]=std::min(val[i],data_[i]);

        }
        return min_val;

    }
    T min_val() {
        T min_val=data_[0];
        for (int i = 0; i < N; ++i) {
            min_val=std::min(min_val,data_[i]);
        }
        return min_val;
    }
    T max_val() {
        T max_val=data_[0];
        for (int i = 0; i < N; ++i) {
            max_val=std::max(max_val,data_[i]);
        }
        return max_val;
    }

}; // end primary template

/******************************************************************************************
 *  Partial Specialization: 2D
 *****************************************************************************************/
template <typename T>
class vec<T, 2> {
public:
    // Named members instead of anonymous struct
    T x, y;
    using value_type = T;
    static constexpr std::size_t dim = 2;
    constexpr vec()  : x(T{}), y(T{}) {}
    constexpr explicit vec(T val)  : x(val), y(val) {}
    constexpr vec(T x_, T y_)  : x(x_), y(y_) {}
    constexpr explicit vec(const std::array<T,2>& arr)  : x(arr[0]), y(arr[1]) {}

    // In-list constructor
    template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) == 2)>>
    constexpr vec(Args... args)  : x{}, y{} {
        T tmp[2] = { (T)args... };
        x = tmp[0]; 
        y = tmp[1];
    }

    /*********************
     *  Operators
     *********************/
    constexpr T& operator[](std::size_t i) {
        switch(i) {
            case 0: return x;
            case 1: return y;
            default: throw std::out_of_range("vec<2> index out of range");
        }
    }
    constexpr const T& operator[](std::size_t i) const {
        switch(i) {
            case 0: return x;
            case 1: return y;
            default: throw std::out_of_range("vec<2> index out of range");
        }
    }

    constexpr bool operator==(const vec<T,2>& other) const  {
        return detail::almost_equal(x, other.x) &&
               detail::almost_equal(y, other.y);
    }
    constexpr bool operator!=(const vec<T,2>& other) const  {
        return !(*this == other);
    }

    constexpr vec<T,2> operator-() const  {
        return vec<T,2>(-x, -y);
    }
    constexpr vec<T,2> operator+(const vec<T,2>& b) const  {
        return {x+b.x, y+b.y};
    }
    constexpr vec<T,2>& operator+=(const vec<T,2>& b)  {
        x += b.x; 
        y += b.y;
        return *this;
    }
    constexpr vec<T,2> operator-(const vec<T,2>& b) const  {
        return {x-b.x, y-b.y};
    }
    constexpr vec<T,2>& operator-=(const vec<T,2>& b)  {
        x -= b.x; 
        y -= b.y;
        return *this;
    }
    constexpr vec<T,2> operator*(T val) const  {
        return {x*val, y*val};
    }
    constexpr vec<T,2>& operator*=(T val)  {
        x*=val; 
        y*=val;
        return *this;
    }
    constexpr vec<T,2> operator/(T val) const {
        if(detail::almost_equal(val, T{})) {
            throw std::invalid_argument("Division by zero");
        }
        return {x/val, y/val};
    }
    constexpr vec<T,2>& operator/=(T val) {
        if(detail::almost_equal(val, T{})) {
            throw std::invalid_argument("Division by zero");
        }
        x/=val; 
        y/=val;
        return *this;
    }
    constexpr vec<T,2> operator*(const vec<T,2>& b) const  {
        return {x*b.x, y*b.y};
    }
    constexpr vec<T,2> operator/(const vec<T,2>& b) const {
        if(detail::almost_equal(b.x, T{}) || detail::almost_equal(b.y, T{})) {
            throw std::invalid_argument("Division by zero in vec<2>");
        }
        return {x/b.x, y/b.y};
    }

    /*********************
     *  Utilities
     *********************/
    constexpr void set(T val)  {
        x = val; y = val;
    }
    constexpr void set(T x_, T y_)  {
        x = x_; y = y_;
    }
    constexpr void set(const std::array<T,2>& arr)  {
        x = arr[0]; y = arr[1];
    }
    constexpr T sqLength() const  {
        return x*x + y*y;
    }
    constexpr T length() const {
        return static_cast<T>(std::sqrt(sqLength()));
    }
    constexpr void unitize() {
        T l = length();
        if(detail::almost_equal(l, T{})) {
            throw std::invalid_argument("Cannot normalize zero-length vector");
        }
        x /= l; 
        y /= l;
    }
    constexpr vec<T,2> unit() const {
        T l = length();
        if(detail::almost_equal(l, T{})) {
            throw std::invalid_argument("Cannot normalize zero-length vector");
        }
        return {x/l, y/l};
    }
    constexpr T dot(const vec<T,2>& other) const  {
        return x*other.x + y*other.y;
    }
    constexpr T distance(const vec<T,2>& other) const {
        return (*this - other).length();
    }
    constexpr T distanceSq(const vec<T,2>& other) const {
        return (*this - other).sqLength();
    }
    constexpr T projection(const vec<T,2>& other) const {
        T dd = other.dot(other);
        if(detail::almost_equal(dd, T{})) {
            throw std::invalid_argument("Projection onto zero-length vector");
        }
        return dot(other)/dd;
    }
    constexpr vec<T,2> project(const vec<T,2>& other) const {
        T dd = other.dot(other);
        if(detail::almost_equal(dd, T{})) {
            throw std::invalid_argument("Projection onto zero-length vector");
        }
        return other*(dot(other)/dd);
    }
    // 2D "cross" => scalar
    constexpr T cross(const vec<T,2>& other) const  {
        return x*other.y - y*other.x;
    }
    constexpr bool collinear(const vec<T,2>& other) const {
        T l1 = length();
        T l2 = other.length();
        if(detail::almost_equal(l1, T{}) || detail::almost_equal(l2, T{})) {
            return false;
        }
        T dotVal = (x/l1)*(other.x/l2) + (y/l1)*(other.y/l2);
        return detail::almost_equal(std::fabs(dotVal), T(1));
    }


    constexpr vec<T,2>max( vec<T,2> &other)  {
        return {std::max(other.x,x),std::max(other.y,y)};

    }
    constexpr vec<T,2>min( vec<T,2> &other)  {
        return {std::min(other.x,x),std::min(other.y,y)};

    }

    friend std::ostream& operator<<(std::ostream& os, const vec<T,2>& v) {
        os << "[" << v.x << "," << v.y << "]";
        return os;
    }
    vec abs() const{

        return {std::abs(x), std::abs(y)};

    }

}; // end partial specialization N=2

/******************************************************************************************
 *  Partial Specialization: 3D
 *****************************************************************************************/
template <typename T>
class vec<T, 3> {
public:
    T x, y, z;
    using value_type = T;
    static constexpr std::size_t dim = 3;
    constexpr vec()  : x(T{}), y(T{}), z(T{}) {}
    constexpr explicit vec(T val)  : x(val), y(val), z(val) {}
    constexpr vec(T x_, T y_, T z_)  : x(x_), y(y_), z(z_) {}
    constexpr explicit vec(const std::array<T,3>& arr)  : x(arr[0]), y(arr[1]), z(arr[2]) {}

    template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) == 3)>>
    constexpr vec(Args... args)  : x{}, y{}, z{} {
        T tmp[3] = { (T)args... };
        x = tmp[0]; 
        y = tmp[1]; 
        z = tmp[2];
    }

    constexpr T& operator[](std::size_t i) {
        switch(i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("vec<3> index out of range");
        }
    }
    constexpr const T& operator[](std::size_t i) const {
        switch(i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("vec<3> index out of range");
        }
    }

    constexpr bool operator==(const vec<T,3>& o) const  {
        return detail::almost_equal(x, o.x) &&
               detail::almost_equal(y, o.y) &&
               detail::almost_equal(z, o.z);
    }
    constexpr bool operator!=(const vec<T,3>& o) const  {
        return !(*this == o);
    }

    constexpr vec<T,3> operator-() const  {
        return {-x, -y, -z};
    }
    constexpr vec<T,3> operator+(const vec<T,3>& b) const  {
        return {x+b.x, y+b.y, z+b.z};
    }
    constexpr vec<T,3>& operator+=(const vec<T,3>& b)  {
        x+=b.x; y+=b.y; z+=b.z; 
        return *this;
    }
    constexpr vec<T,3> operator-(const vec<T,3>& b) const  {
        return {x-b.x, y-b.y, z-b.z};
    }
    constexpr vec<T,3>& operator-=(const vec<T,3>& b)  {
        x-=b.x; y-=b.y; z-=b.z;
        return *this;
    }
    constexpr vec<T,3> operator*(T val) const  {
        return {x*val, y*val, z*val};
    }
    constexpr vec<T,3>& operator*=(T val)  {
        x*=val; y*=val; z*=val;
        return *this;
    }
    constexpr vec<T,3> operator/(T val) const {
        if(detail::almost_equal(val, T{})) {
            throw std::invalid_argument("Division by zero");
        }
        return {x/val, y/val, z/val};
    }
    constexpr vec<T,3>& operator/=(T val) {
        if(detail::almost_equal(val, T{})) {
            throw std::invalid_argument("Division by zero");
        }
        x/=val; y/=val; z/=val;
        return *this;
    }
    constexpr vec<T,3> operator*(const vec<T,3>& b) const  {
        return {x*b.x, y*b.y, z*b.z};
    }
    constexpr vec<T,3> operator/(const vec<T,3>& b) const {
        if(detail::almost_equal(b.x, T{}) ||
           detail::almost_equal(b.y, T{}) ||
           detail::almost_equal(b.z, T{})) {
            throw std::invalid_argument("Division by zero in vec<3>");
        }
        return {x/b.x, y/b.y, z/b.z};
    }

    // Utilities
    constexpr void set(T val)  { x=val; y=val; z=val; }
    constexpr void set(T xx, T yy, T zz)  { x=xx; y=yy; z=zz; }
    constexpr void set(const std::array<T,3>& arr)  {
        x=arr[0]; y=arr[1]; z=arr[2];
    }

    constexpr T sqLength() const  { return x*x + y*y + z*z; }
    constexpr T length() const { return static_cast<T>(std::sqrt(sqLength())); }
    constexpr void unitize() {
        T l = length();
        if(detail::almost_equal(l, T{})) {
            throw std::invalid_argument("Cannot normalize zero-length vector");
        }
        x/=l; y/=l; z/=l;
    }
    constexpr vec<T,3> unit() const {
        T l = length();
        if(detail::almost_equal(l, T{})) {
            throw std::invalid_argument("Cannot normalize zero-length vector");
        }
        return {x/l, y/l, z/l};
    }
    constexpr T dot(const vec<T,3>& b) const  {
        return x*b.x + y*b.y + z*b.z;
    }
    constexpr vec<T,3> cross(const vec<T,3>& b) const  {
        return {
            y*b.z - z*b.y,
            z*b.x - x*b.z,
            x*b.y - y*b.x
        };
    }
    constexpr T distance(const vec<T,3>& other) const {
        return (*this - other).length();
    }
    constexpr T distanceSq(const vec<T,3>& other) const  {
        return (*this - other).sqLength();
    }
    constexpr T projection(const vec<T,3>& other) const {
        T dd = other.dot(other);
        if(detail::almost_equal(dd, T{})) {
            throw std::invalid_argument("Projection onto zero-length vector");
        }
        return dot(other)/dd;
    }
    constexpr vec<T,3> project(const vec<T,3>& other) const {
        T dd = other.dot(other);
        if(detail::almost_equal(dd, T{})) {
            throw std::invalid_argument("Projection onto zero-length vector");
        }
        return other*(dot(other)/dd);
    }
    constexpr bool collinear(const vec<T,3>& other) const {
        T l1 = length();
        T l2 = other.length();
        if(detail::almost_equal(l1, T{}) || detail::almost_equal(l2, T{})) {
            return false;
        }
        T dotVal = (x/l1)*(other.x/l2) + (y/l1)*(other.y/l2) + (z/l1)*(other.z/l2);
        return detail::almost_equal(std::fabs(dotVal), T(1));
    }
    constexpr vec<T,3> max( vec<T,3> &other)  {
        return {std::max(other.x,x),std::max(other.y,y),std::max(other.z,z)};

    }
    constexpr vec<T,3> min( vec<T,3> &other)  {
        return {std::min(other.x,x),std::min(other.y,y),std::min(other.z,z)};

    }
    friend std::ostream& operator<<(std::ostream& os, const vec<T,3>& v) {
        os << "[" << v.x << "," << v.y << "," << v.z << "]";
        return os;
    }
    vec abs() const{

        return {std::abs(x), std::abs(y), std::abs(z)};

    }

    bool zero()const {
        return (x==0)&y==0&z==0;
    }

};

/******************************************************************************************
 *  Partial Specialization: 4D
 *****************************************************************************************/
template <typename T>
class vec<T, 4> {
public:
    T x, y, z, w;
    static constexpr std::size_t dim = 4;
    using value_type = T;
    /*********************
     *  Constructors
     *********************/
    constexpr vec()  : x(T{}), y(T{}), z(T{}), w(T{}) {}
    constexpr explicit vec(T val)  : x(val), y(val), z(val), w(val) {}
    constexpr vec(T xx, T yy, T zz, T ww)  : x(xx), y(yy), z(zz), w(ww) {}
    constexpr explicit vec(const std::array<T,4>& arr) 
        : x(arr[0]), y(arr[1]), z(arr[2]), w(arr[3]) {}

    template <typename... Args, typename = std::enable_if_t<(sizeof...(Args) == 4)>>
    constexpr vec(Args... args)  : x{}, y{}, z{}, w{} {
        T tmp[4] = { (T)args... };
        x = tmp[0]; 
        y = tmp[1]; 
        z = tmp[2]; 
        w = tmp[3];
    }

    constexpr T& operator[](std::size_t i) {
        switch(i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: throw std::out_of_range("vec<4> index out of range");
        }
    }
    constexpr const T& operator[](std::size_t i) const {
        switch(i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: throw std::out_of_range("vec<4> index out of range");
        }
    }

    constexpr bool operator==(const vec<T,4>& o) const  {
        return detail::almost_equal(x, o.x) &&
               detail::almost_equal(y, o.y) &&
               detail::almost_equal(z, o.z) &&
               detail::almost_equal(w, o.w);
    }
    constexpr bool operator!=(const vec<T,4>& o) const  {
        return !(*this == o);
    }

    constexpr vec<T,4> operator-() const  {
        return {-x, -y, -z, -w};
    }
    constexpr vec<T,4> operator+(const vec<T,4>& b) const  {
        return {x+b.x, y+b.y, z+b.z, w+b.w};
    }
    constexpr vec<T,4>& operator+=(const vec<T,4>& b)  {
        x+=b.x; y+=b.y; z+=b.z; w+=b.w;
        return *this;
    }
    constexpr vec<T,4> operator-(const vec<T,4>& b) const  {
        return {x-b.x, y-b.y, z-b.z, w-b.w};
    }
    constexpr vec<T,4>& operator-=(const vec<T,4>& b)  {
        x-=b.x; y-=b.y; z-=b.z; w-=b.w;
        return *this;
    }
    constexpr vec<T,4> operator*(T val) const  {
        return {x*val, y*val, z*val, w*val};
    }
    constexpr vec<T,4>& operator*=(T val)  {
        x*=val; y*=val; z*=val; w*=val;
        return *this;
    }
    constexpr vec<T,4> operator/(T val) const {
        if(detail::almost_equal(val, T{})) {
            throw std::invalid_argument("Division by zero");
        }
        return {x/val, y/val, z/val, w/val};
    }
    constexpr vec<T,4>& operator/=(T val) {
        if(detail::almost_equal(val, T{})) {
            throw std::invalid_argument("Division by zero");
        }
        x/=val; y/=val; z/=val; w/=val;
        return *this;
    }

    /// Component-wise multiply
    constexpr vec<T,4> operator*(const vec<T,4>& b) const  {
        return {x*b.x, y*b.y, z*b.z, w*b.w};
    }
    /// Component-wise divide
    constexpr vec<T,4> operator/(const vec<T,4>& b) const {
        if(detail::almost_equal(b.x, T{}) ||
           detail::almost_equal(b.y, T{}) ||
           detail::almost_equal(b.z, T{}) ||
           detail::almost_equal(b.w, T{})) {
            throw std::invalid_argument("Division by zero in vec<4>");
        }
        return {x/b.x, y/b.y, z/b.z, w/b.w};
    }

    /*********************
     *  Utilities
     *********************/
    constexpr void set(T val)  {
        x=val; y=val; z=val; w=val;
    }
    constexpr void set(T xx, T yy, T zz, T ww)  {
        x=xx; y=yy; z=zz; w=ww;
    }
    constexpr void set(const std::array<T,4>& arr)  {
        x=arr[0]; y=arr[1]; z=arr[2]; w=arr[3];
    }
    constexpr T sqLength() const  {
        return x*x + y*y + z*z + w*w;
    }
    constexpr T length() const {
        return static_cast<T>(std::sqrt(sqLength()));
    }

    /// "Homogeneous" unitize approach:
    /// some treatments treat w as 1 always. We keep consistent with the original code:
    /// x/(w*l), y/(w*l), z/(w*l), w=1
    constexpr void unitize() {
        T l = length();
        if(detail::almost_equal(l, T{})) {
            throw std::invalid_argument("Cannot normalize zero-length vector (4D)");
        }
        // Guard against zero w so as not to multiply by 1/w -> INF
        if(detail::almost_equal(w, T{})) {
            throw std::invalid_argument("Cannot normalize 4D vector with w=0 in this design");
        }
        x /= (w*l);
        y /= (w*l);
        z /= (w*l);
        w  = T(1);
    }
    constexpr vec<T,4> unit() const {
        vec<T,4> result(*this);
        result.unitize();
        return result;
    }

    /// Convert to a 3D vector by dividing each x,y,z by w
    constexpr vec<T,3> to_vec3() const {
        if(detail::almost_equal(w, T{})) {
            throw std::invalid_argument("Cannot convert 4D to 3D with w=0");
        }
        return vec<T,3>(x/w, y/w, z/w);
    }

    /// Dot product ignoring w in original code or not?
    /// Original code did: return x*b.x + y*b.y + z*b.z ( it omitted w!? )
    /// Let's fix it to be consistent with a 4D dot:
    constexpr T dot(const vec<T,4>& b) const  {
        // We'll do the full 4D dot product
        return x*b.x + y*b.y + z*b.z + w*b.w;
    }

    /// Cross product in 4D is not standard, but the original code tries
    /// to do it by converting to vec3. We fix the bug:
    ///   result.x = v1.y * v2.z - v1.z * v2.y, ...
    /// This is purely a 3D cross ignoring w, returning a 4D vector with w=1.
    constexpr vec<T,4> cross(const vec<T,4>& b) const {
        // interpret x/w, y/w, z/w
        if(detail::almost_equal(w, T{}) || detail::almost_equal(b.w, T{})) {
            throw std::invalid_argument("Cannot compute cross for 4D vectors with w=0");
        }
        vec<T,3> v1(x/w, y/w, z/w);
        vec<T,3> v2(b.x/b.w, b.y/b.w, b.z/b.w);
        auto c = v1.cross(v2);
        return vec<T,4>(c.x, c.y, c.z, T(1));
    }

    /// Distance to another 4D vector
    constexpr T distance(const vec<T,4>& other) const {
        return (*this - other).length();
    }
    constexpr T distanceSq(const vec<T,4>& other) const {
        return (*this - other).sqLength();
    }

    /// Projection length: (this·other)/(other·other)
    constexpr T projection(const vec<T,4>& other) const {
        T dd = other.dot(other);
        if(detail::almost_equal(dd, T{})) {
            throw std::invalid_argument("Projection onto zero-length 4D vector");
        }
        return dot(other)/dd;
    }

    /// vec projection: other * ((this·other)/(other·other))
    constexpr vec<T,4> project(const vec<T,4>& other) const {
        T dd = other.dot(other);
        if(detail::almost_equal(dd, T{})) {
            throw std::invalid_argument("Projection onto zero-length 4D vector");
        }
        return other*(dot(other)/dd);
    }
    constexpr bool collinear(const vec<T,4>& other) const {
        // check if 3D parts are collinear once we interpret them
        // as x/w, y/w, z/w. We do the same for 'other'.
        if(detail::almost_equal(w, T{}) || detail::almost_equal(other.w, T{})) {
            return false;
        }
        vec<T,3> v1(x/w, y/w, z/w);
        vec<T,3> v2(other.x/other.w, other.y/other.w, other.z/other.w);
        return v1.collinear(v2);
    }

    friend std::ostream& operator<<(std::ostream& os, const vec<T,4>& v) {
        os << "[" << v.x << "," << v.y << "," << v.z << "," << v.w << "]";
        return os;
    }

    vec abs() const{

        return {std::abs(x), std::abs(y), std::abs(z), std::abs(w)};

    }
}; // end partial specialization N=4

/******************************************************************************************
 *  Free Function Operators (Scalar * vec)
 *****************************************************************************************/
template <typename T, std::size_t N>
constexpr vec<T,N> operator*(T scalar, const vec<T,N>& vec)  {
    return vec * scalar;
}

/******************************************************************************************
 *  Cartesian <-> Spherical specialized for double 2D & 3D
 *****************************************************************************************/

// cartesian_to_spherical(3D -> 3D): 
//   r = sqrt(x^2 + y^2 + z^2)
//   θ = atan2( sqrt(x^2 + y^2), z )
//   φ = atan2(y, x)
inline void cartesian_to_spherical(const vec<double,3>& xyz, vec<double,3>& rtp) {
    double x2y2 = xyz.x*xyz.x + xyz.y*xyz.y;
    rtp[0] = std::sqrt(x2y2 + xyz.z*xyz.z);       // r
    rtp[1] = std::atan2(std::sqrt(x2y2), xyz.z);  // θ
    rtp[2] = std::atan2(xyz.y, xyz.x);            // φ
}

// spherical_to_cartesian(3D -> 3D):
//   x = r * sin(θ) * cos(φ)
//   y = r * sin(θ) * sin(φ)
//   z = r * cos(θ)
inline void spherical_to_cartesian(const vec<double,3>& rtp, vec<double,3>& xyz) {
    double r = rtp[0];
    double theta = rtp[1];
    double phi = rtp[2];
    xyz[0] = r*std::sin(theta)*std::cos(phi);
    xyz[1] = r*std::sin(theta)*std::sin(phi);
    xyz[2] = r*std::cos(theta);
}

// cartesian_to_spherical(3D -> 2D): 
//   θ=atan2( sqrt(x^2 + y^2), z )
//   φ=atan2(y, x)
inline void cartesian_to_spherical(const vec<double,3>& xyz, vec<double,2>& tp) {
    double x2y2 = xyz.x*xyz.x + xyz.y*xyz.y;
    tp[0] = std::atan2(std::sqrt(x2y2), xyz.z);  // θ
    tp[1] = std::atan2(xyz.y, xyz.x);            // φ
}

// spherical_to_cartesian(2D -> 3D) assuming r=1
inline void spherical_to_cartesian(const vec<double,2>& tp, vec<double,3>& xyz) {
    double r = 1.0;
    double theta = tp[0];
    double phi   = tp[1];
    xyz[0] = r*std::sin(theta)*std::cos(phi);
    xyz[1] = r*std::sin(theta)*std::sin(phi);
    xyz[2] = r*std::cos(theta);
}

// cartesian_to_spherical over std::vector
inline void cartesian_to_spherical(const std::vector<vec<double,3>>& xyz,
                                   std::vector<vec<double,2>>& tp) {
    tp.resize(xyz.size());
    for(std::size_t i = 0; i < xyz.size(); ++i) {
        double x2y2 = xyz[i].x*xyz[i].x + xyz[i].y*xyz[i].y;
        tp[i][0] = std::atan2(std::sqrt(x2y2), xyz[i].z);
        tp[i][1] = std::atan2(xyz[i].y, xyz[i].x);
    }
}
template<typename T,int N>
int max_axis(const vec<T,N> &v) {
    T max_val =std::numeric_limits<T>::min();
    size_t max_idx=0;
    for(size_t i=0;i<N;++i) {

    if (v[i] > max_val) {
        max_val=v[i];
        max_idx=i;
    }
    }

    return max_idx;
    };


    template<typename T>
    inline int max_axis(const vec<T,2> &v) {
    // For a 3D vector v = (v[0], v[1], v[2])
    // return index of the largest component
    if (v[0] > v[1]) return 0;
    return 1;

}

    template<typename T>
    inline int max_axis(const vec<T,3> &v) {
    // For a 3D vector v = (v[0], v[1], v[2])
    // return index of the largest component
    if (v[0] > v[1] && v[0] > v[2]) return 0;
    if (v[1] > v[2]) return 1;
    return 2;
}



/******************************************************************************************
 *  End of single-header vector library
 *****************************************************************************************/

    using vec2=vec<double,2>;
    using vec3=vec<double,3>;
    using vec4=vec<double,4>;

    using vec2d=vec<double,2>;
    using vec3d=vec<double,3>;
    using vec4d=vec<double,4>;

    using vec2f=vec<float,2>;
    using vec3f=vec<float,3>;
    using vec4f=vec<float,4>;
} // end namespace bvh

#endif // VEC_H
