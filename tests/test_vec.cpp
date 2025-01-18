//
// Created by Andrew Astakhov on 14.01.25.
//
/*****************************************************************************************
 *  Simple Test Suite for the Refactored vec Library
 *  (No external test framework; just basic asserts and try/catch blocks)
 *
 *  Usage:
 *      1. Make sure you have the refactored "Vec.h" in the same directory.
 *      2. Compile with a C++ compiler, e.g.:
 *         g++ -std=c++17 -o test_vec test_vec.cpp
 *      3. Run: ./test_vec
 *
 *  Note: This test suite is not exhaustive, but it covers the major functionalities
 *        including arithmetic, dot/cross, lengths, projections, and exception cases.
 *****************************************************************************************/

#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include "../src/vec.h"  // Include your refactored single-header vector library
#include "src/aabb.h"
#include "src/bvh.h"
using namespace bvh;

// Helper: floating-point approximate comparison for tests.
inline bool almost_equal(double a, double b, double eps = 1e-12)
{
    return std::fabs(a - b) < eps;
}

/*****************************************************************************************
 *  Tests for 2D vectors
 *****************************************************************************************/
static void testvec2D()
{
    // Constructors
    vec<double,2> v1;                // default ctor
    assert(almost_equal(v1.x, 0.0) && almost_equal(v1.y, 0.0));

    vec<double,2> v2(5.0);           // fill constructor
    assert(almost_equal(v2.x, 5.0) && almost_equal(v2.y, 5.0));

    vec<double,2> v3(1.0, 2.0);      // direct constructor
    assert(almost_equal(v3.x, 1.0) && almost_equal(v3.y, 2.0));

    std::array<double,2> arr2 = {3.0, 4.0};
    vec<double,2> v4(arr2);          // from std::array
    assert(almost_equal(v4.x, 3.0) && almost_equal(v4.y, 4.0));

    // Basic arithmetic
    auto sum = v3 + v4;                 // (1,2) + (3,4) = (4,6)
    assert(almost_equal(sum.x, 4.0) && almost_equal(sum.y, 6.0));

    auto diff = v4 - v3;                // (3,4) - (1,2) = (2,2)
    assert(almost_equal(diff.x, 2.0) && almost_equal(diff.y, 2.0));

    auto neg = -v3;                     // -(1,2) = (-1,-2)
    assert(almost_equal(neg.x, -1.0) && almost_equal(neg.y, -2.0));

    // Scalar multiply/divide
    auto scaled = v3 * 2.0;             // (1,2)*2 = (2,4)
    assert(almost_equal(scaled.x, 2.0) && almost_equal(scaled.y, 4.0));

    scaled /= 2.0;                      // (2,4)/2 = (1,2)
    assert(almost_equal(scaled.x, 1.0) && almost_equal(scaled.y, 2.0));

    // Component-wise multiply/divide
    vec<double,2> cwMul = v3 * v4;   // (1,2)*(3,4) = (3,8)
    assert(almost_equal(cwMul.x, 3.0) && almost_equal(cwMul.y, 8.0));

    vec<double,2> cwDiv = v4 / v3;   // (3,4)/(1,2) = (3,2)
    assert(almost_equal(cwDiv.x, 3.0) && almost_equal(cwDiv.y, 2.0));

    // Dot, length, distance
    double dotVal = v3.dot(v4);         // (1,2)·(3,4) = 1*3 + 2*4 = 11
    assert(almost_equal(dotVal, 11.0));

    assert(almost_equal(v3.sqLength(), 5.0));  // (1,2).sqLength() = 1 + 4 = 5
    assert(almost_equal(v3.length(), std::sqrt(5.0)));

    double dist = v3.distance(v4);      // distance((1,2),(3,4))=sqrt(2^2+2^2)=sqrt(8)=2.828...
    assert(almost_equal(dist, std::sqrt(8.0)));

    // Cross in 2D is a scalar
    auto crossVal = v3.cross(v4);       // cross((1,2),(3,4))=1*4 - 2*3=4 -6= -2
    assert(almost_equal(crossVal, -2.0));

    // Projection
    auto projLength = v3.projection(v4);  // (1,2)·(3,4)/ (3,4)·(3,4)=11/25=0.44
    assert(almost_equal(projLength, 11.0/25.0));

    auto projVec = v3.project(v4);        // = v4 * 0.44 approx => (3*0.44,4*0.44)...
    assert(almost_equal(projVec.x, 3.0 * 0.44));
    assert(almost_equal(projVec.y, 4.0 * 0.44));

    // Collinearity
    vec<double,2> vCol1(2.0, 4.0);
    assert(vCol1.collinear(vec<double,2>(1.0,2.0)) == true);

    // Test normalization
    vec<double,2> vNorm(3.0, 4.0);
    vNorm.unitize();
    assert(almost_equal(vNorm.length(), 1.0));
    assert(almost_equal(vNorm.x, 3.0/5.0) && almost_equal(vNorm.y, 4.0/5.0));

    // Test exceptions
    {
        bool thrown = false;
        try {
            vec<double,2> zeroVec(0.0,0.0);
            zeroVec.unitize();  // should throw
        } catch(const std::invalid_argument&) {
            thrown = true;
        }
        assert(thrown && "Expected exception for normalizing zero vector in 2D");
    }

    {
        bool thrown = false;
        try {
            vec<double,2> testDiv(1.0,2.0);
            testDiv /= 0.0;  // should throw
        } catch(const std::invalid_argument&) {
            thrown = true;
        }
        assert(thrown && "Expected exception for dividing by zero in 2D");
    }

    std::cout << "testvec2D() passed.\n";
}

/*****************************************************************************************
 *  Tests for 3D vectors
 *****************************************************************************************/
static void testvec3D()
{
    // Constructors
    vec<double,3> v1;               // (0,0,0)
    assert(almost_equal(v1.x,0.0)&&almost_equal(v1.y,0.0)&&almost_equal(v1.z,0.0));

    vec<double,3> v2(2.0);          // (2,2,2)
    assert(almost_equal(v2.x,2.0)&&almost_equal(v2.y,2.0)&&almost_equal(v2.z,2.0));

    vec<double,3> v3(1.0, 2.0, 3.0);
    assert(almost_equal(v3.x,1.0)&&almost_equal(v3.y,2.0)&&almost_equal(v3.z,3.0));

    // Arithmetic
    auto add = v2 + v3;                // (2,2,2)+(1,2,3)=(3,4,5)
    assert(almost_equal(add.x,3.0)&&almost_equal(add.y,4.0)&&almost_equal(add.z,5.0));

    auto sub = v3 - v2;                // (1,2,3)-(2,2,2)=( -1,0,1 )
    assert(almost_equal(sub.x,-1.0)&&almost_equal(sub.y,0.0)&&almost_equal(sub.z,1.0));

    auto neg = -v3;                    // -(1,2,3)=(-1,-2,-3)
    assert(almost_equal(neg.x,-1.0)&&almost_equal(neg.y,-2.0)&&almost_equal(neg.z,-3.0));

    auto scaled = v3 * 3.0;            // (1,2,3)*3 = (3,6,9)
    assert(almost_equal(scaled.x,3.0)&&almost_equal(scaled.y,6.0)&&almost_equal(scaled.z,9.0));

    scaled /= 3.0;                     // (3,6,9)/3 = (1,2,3)
    assert(almost_equal(scaled.x,1.0)&&almost_equal(scaled.y,2.0)&&almost_equal(scaled.z,3.0));

    // Dot & Cross
    vec<double,3> vA(1,0,0), vB(0,1,0);
    double dotVal = vA.dot(vB);        // (1,0,0)·(0,1,0)=0
    assert(almost_equal(dotVal,0.0));

    auto crossVal = vA.cross(vB);      // cross((1,0,0),(0,1,0))=(0,0,1)
    assert(almost_equal(crossVal.x,0.0)&&almost_equal(crossVal.y,0.0)&&almost_equal(crossVal.z,1.0));

    // Length / Distance
    vec<double,3> vC(3.0,4.0,12.0);
    assert(almost_equal(vC.sqLength(), 3.0*3.0+4.0*4.0+12.0*12.0));
    assert(almost_equal(vC.length(), 13.0));

    auto dist = v3.distance(vC); // distance((1,2,3),(3,4,12))= sqrt( (3-1)^2+(4-2)^2+(12-3)^2 )
    assert(almost_equal(dist, std::sqrt(2*2 + 2*2 + 9*9)));  // sqrt(4+4+81)= sqrt(89)

    // Projection
    vec<double,3> vD(3.0,0.0,0.0);
    auto projLength = v3.projection(vD);  // = dot((1,2,3),(3,0,0)) / dot((3,0,0),(3,0,0)) = 3/9=1/3
    assert(almost_equal(projLength, 1.0/3.0));

    // Collinearity
    vec<double,3> vCol1(2.0,4.0,6.0);
    assert(v3.collinear(vCol1) == true);  // (1,2,3) vs (2,4,6)

    // Normalization
    vec<double,3> vNorm(1.0,2.0,2.0);
    vNorm.unitize();
    assert(almost_equal(vNorm.length(),1.0));

    // Exceptions
    {
        bool thrown = false;
        try {
            vec<double,3> zeroVec(0.0,0.0,0.0);
            zeroVec.unitize();
        } catch(const std::invalid_argument&) {
            thrown = true;
        }
        assert(thrown && "Expected exception for normalizing zero vector in 3D");
    }

    {
        bool thrown = false;
        try {
            vec<double,3> testDiv(1.0,2.0,3.0);
            testDiv /= 0.0;
        } catch(const std::invalid_argument&) {
            thrown = true;
        }
        assert(thrown && "Expected exception for dividing by zero in 3D");
    }

    std::cout << "testvec3D() passed.\n";
}

/*****************************************************************************************
 *  Tests for 4D vectors
 *****************************************************************************************/
static void testvec4D()
{
    // Constructors
    vec<double,4> v1; // (0,0,0,0)
    assert(almost_equal(v1.x,0.0)&&almost_equal(v1.y,0.0)&&almost_equal(v1.z,0.0)&&almost_equal(v1.w,0.0));

    vec<double,4> v2(1.0); // (1,1,1,1)
    assert(almost_equal(v2.x,1.0)&&almost_equal(v2.y,1.0)&&almost_equal(v2.z,1.0)&&almost_equal(v2.w,1.0));

    vec<double,4> v3(1.0,2.0,3.0,4.0);
    assert(almost_equal(v3.x,1.0)&&almost_equal(v3.y,2.0)&&almost_equal(v3.z,3.0)&&almost_equal(v3.w,4.0));

    // Arithmetic
    auto add = v2 + v3; // (1,1,1,1)+(1,2,3,4)=(2,3,4,5)
    assert(almost_equal(add.x,2.0)&&almost_equal(add.y,3.0)&&almost_equal(add.z,4.0)&&almost_equal(add.w,5.0));

    auto sub = v3 - v2; // (1,2,3,4)-(1,1,1,1)=(0,1,2,3)
    assert(almost_equal(sub.x,0.0)&&almost_equal(sub.y,1.0)&&almost_equal(sub.z,2.0)&&almost_equal(sub.w,3.0));

    auto neg = -v3;     // -(1,2,3,4)=(-1,-2,-3,-4)
    assert(almost_equal(neg.x,-1.0)&&almost_equal(neg.y,-2.0)&&almost_equal(neg.z,-3.0)&&almost_equal(neg.w,-4.0));

    auto scaled = v3 * 2.0; // (1,2,3,4)*2=(2,4,6,8)
    assert(almost_equal(scaled.x,2.0)&&almost_equal(scaled.y,4.0)&&almost_equal(scaled.z,6.0)&&almost_equal(scaled.w,8.0));

    scaled /= 2.0; // (2,4,6,8)/2=(1,2,3,4)
    assert(almost_equal(scaled.x,1.0)&&almost_equal(scaled.y,2.0)&&almost_equal(scaled.z,3.0)&&almost_equal(scaled.w,4.0));

    // Dot product
    double dotVal = v2.dot(v3); // (1,1,1,1)·(1,2,3,4)=1+2+3+4=10
    assert(almost_equal(dotVal,10.0));

    // Convert to vec3
    {
        bool thrown = false;
        try {
            auto as3 = v3.to_vec3(); // (1/4,2/4,3/4)
            assert(almost_equal(as3.x, 1.0/4.0));
            assert(almost_equal(as3.y, 2.0/4.0));
            assert(almost_equal(as3.z, 3.0/4.0));
        } catch(...) {
            // Should NOT throw
            thrown = true;
        }
        assert(!thrown && "Unexpected exception in to_vec3(4D->3D)");
    }

    // Cross product in 4D (via 3D part) 
    // e.g. (1,0,0,1) x (0,1,0,1) => interpret as 3D(1,0,0) x 3D(0,1,0) => (0,0,1), w=1
    {
        vec<double,4> vx(1,0,0,1);
        vec<double,4> vy(0,1,0,1);
        auto crossVal = vx.cross(vy);
        assert(almost_equal(crossVal.x,0.0));
        assert(almost_equal(crossVal.y,0.0));
        assert(almost_equal(crossVal.z,1.0));
        assert(almost_equal(crossVal.w,1.0));
    }

    // Unitize: (1,2,3,4), length= sqrt(1+4+9+16)= sqrt(30)=5.4772..., w=4
    //   -> after unitize => x/(4*l), y/(4*l), z/(4*l), w=1
    {
        vec<double,4> temp(1,2,3,4);
        double l = temp.length(); // sqrt(30)
        temp.unitize();
        // x'=1/(4*sqrt(30)), y'=2/(4*sqrt(30)), z'=3/(4*sqrt(30)), w'=1
        assert(almost_equal(temp.x, 1.0/(4.0*l)));
        assert(almost_equal(temp.y, 2.0/(4.0*l)));
        assert(almost_equal(temp.z, 3.0/(4.0*l)));
        assert(almost_equal(temp.w, 1.0));
    }

    // Projection
    {
        vec<double,4> a(1,2,0,1);
        vec<double,4> b(2,0,2,1);
        double projLen = a.projection(b);
        // dot(a,b) = [ (1,2,0,1)·(2,0,2,1) ] = 1*2 + 2*0 + 0*2 + 1*1=2+1=3
        // dot(b,b) = (2,0,2,1)·(2,0,2,1)=4+4+1=9
        // => 3/9=1/3=0.3333
        assert(almost_equal(projLen, 1.0/3.0));

        auto projVec = a.project(b);
        // => b*(1/3)= (2,0,2,1)*(1/3)= (2/3,0,2/3,1/3)
        assert(almost_equal(projVec.x, 2.0/3.0));
        assert(almost_equal(projVec.y, 0.0));
        assert(almost_equal(projVec.z, 2.0/3.0));
        assert(almost_equal(projVec.w, 1.0/3.0));
    }

    // Exceptions
    {
        bool thrown = false;
        try {
            vec<double,4> zeroW(1,2,3,0);
            zeroW.unitize(); // w=0 => throw
        } catch(const std::invalid_argument&) {
            thrown = true;
        }
        assert(thrown && "Expected exception for w=0 in 4D unitize");
    }

    {
        bool thrown = false;
        try {
            vec<double,4> allZero(0,0,0,0);
            allZero.unitize(); // length=0 => throw
        } catch(const std::invalid_argument&) {
            thrown = true;
        }
        assert(thrown && "Expected exception for zero-length in 4D");
    }

    std::cout << "testvec4D() passed.\n";
}

/*****************************************************************************************
 *  Tests for Cartesian <-> Spherical Conversions
 *****************************************************************************************/
static void testCartesianSpherical()
{
    // 3D -> spherical -> 3D
    {
        vec<double,3> cart(3.0, 4.0, 12.0);
        vec<double,3> sph;
        cartesian_to_spherical(cart, sph);
        // r=13, θ=atan2(sqrt(3^2+4^2), 12)=atan2(5,12), φ=atan2(4,3)
        assert(almost_equal(sph[0], 13.0));
        assert(almost_equal(std::tan(sph[1]), 5.0/12.0));
        assert(almost_equal(std::tan(sph[2]), 4.0/3.0));

        vec<double,3> cart2;
        spherical_to_cartesian(sph, cart2);
        // should be back to (3,4,12) (within floating tolerance)
        assert(almost_equal(cart2.x,3.0));
        assert(almost_equal(cart2.y,4.0));
        assert(almost_equal(cart2.z,12.0));
    }

    // 3D -> 2D
    {
        vec<double,3> cart(1.0, 1.0, 0.0);
        vec<double,2> sph2;
        cartesian_to_spherical(cart, sph2);
        // => θ=atan2(sqrt(1+1),0)=atan2(sqrt(2),0)= pi/2
        //    φ=atan2(1,1)= pi/4
        assert(almost_equal(sph2[0], M_PI/2));
        assert(almost_equal(sph2[1], M_PI/4));
    }

    // 2D -> 3D (assuming r=1)
    {
        vec<double,2> angles(M_PI/2, 0);
        vec<double,3> out;
        spherical_to_cartesian(angles, out);
        // => (sin(pi/2)*cos(0), sin(pi/2)*sin(0), cos(pi/2)) = (1,0,0)
        assert(almost_equal(out.x,1.0));
        assert(almost_equal(out.y,0.0));
        assert(almost_equal(out.z,0.0));
    }

    // Batch cart->spherical
    {
        std::vector<vec<double,3>> inVecs({
            {3.0, 4.0, 12.0},
            {0.0, 1.0, 0.0},
            {10.0, 0.0, 10.0}
        });
        std::vector<vec<double,2>> outVecs;
        cartesian_to_spherical(inVecs, outVecs);
        assert(outVecs.size()==3);
        // Just do a quick check on one:
        // second element => (0,1,0) => θ=atan2(sqrt(0+1),0)=atan2(1,0)= pi/2
        //                              φ=atan2(1,0)= pi/2
        assert(almost_equal(outVecs[1][0], M_PI/2));
        assert(almost_equal(outVecs[1][1], M_PI/2));
    }

    std::cout << "testCartesianSpherical() passed.\n";
}

/*****************************************************************************************
 *  Main test runner
 *****************************************************************************************/
int main()
{
    try {
        testvec2D();
        testvec3D();
        testvec4D();
        testCartesianSpherical();
    } 
    catch(const std::exception& e) {
        std::cerr << "Test suite failed with exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "All tests passed successfully!\n";
    return 0;
}
