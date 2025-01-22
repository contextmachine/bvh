//
// Created by Andrew Astakhov on 14.01.25.
//
/*****************************************************************************************
 *  Simple Test Suite for the AABB Class
 *  (No external test framework; just basic asserts and try/catch blocks)
 *
 *  Usage:
 *      1. Make sure you have the "aabb.h" in the same directory.
 *      2. Compile with a C++ compiler, e.g.:
 *         g++ -std=c++17 -o test_aabb test_aabb.cpp
 *      3. Run: ./test_aabb
 *
 *  Note: This test suite is not exhaustive, but it covers the major functionalities
 *        including construction, intersection, volume, and expansion.
 *****************************************************************************************/

#include <iostream>
#include <cassert>
#include "bvh/vec.h"
#include "bvh/aabb.h"
using namespace bvh;

// Helper: floating-point approximate comparison for tests.
inline bool almost_equal(double a, double b, double eps = 1e-12)
{
    return std::fabs(a - b) < eps;
}

/*****************************************************************************************
 *  Tests for AABB
 *****************************************************************************************/
static void testAABB()
{
    // Test default constructor
    AABB<vec<double, 3>> aabb1;
    assert(almost_equal(aabb1.min.x, std::numeric_limits<double>::max()));
    assert(almost_equal(aabb1.max.x, std::numeric_limits<double>::lowest()));

    // Test constructor with min and max
    vec<double, 3> minVec(0.0, 0.0, 0.0);
    vec<double, 3> maxVec(1.0, 1.0, 1.0);
    AABB aabb2(minVec, maxVec);
    assert(almost_equal(aabb2.min.x, 0.0));
    assert(almost_equal(aabb2.max.x, 1.0));

    // Test intersection
    AABB aabb3(vec<double, 3>(0.5, 0.5, 0.5), vec<double, 3>(1.5, 1.5, 1.5));
    AABB<vec<double,3>> result;
    assert(aabb2.intersection(aabb3, result));
    assert(almost_equal(result.min.x, 0.5));
    assert(almost_equal(result.max.x, 1.0));

    // Test volume
    assert(almost_equal(aabb2.volume(), 1.0));

    // Test expand by point
    aabb2.expand(vec<double, 3>(2.0, 2.0, 2.0));
    assert(almost_equal(aabb2.max.x, 2.0));

    std::cout << "testAABB() passed.\n";
}

/*****************************************************************************************
 *  Main test runner
 *****************************************************************************************/
int main()
{
    try {
        testAABB();
    } 
    catch(const std::exception& e) {
        std::cerr << "Test suite failed with exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "All tests passed successfully!\n";
    return 0;
} 