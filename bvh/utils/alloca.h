//
// Created by Andrew Astakhov on 28.01.25.
//

#ifndef BVH_UTILS_ALLOCA_H
#define BVH_UTILS_ALLOCA_H

#include <stdio.h>

// Use _alloca only for MSVC on Windows
#if defined(_MSC_VER)
    #include <malloc.h> // Required for alloca/_alloca
    #define alloca _alloca
#endif

#endif //BVH_UTILS_ALLOCA_H
