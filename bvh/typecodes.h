//
// Created by Andrew Astakhov on 23.01.25.
//

#ifndef TYPECODES_H
#define TYPECODES_H


#include <cstdint>

namespace bvh {
    /**********************************************************************
     *  Helper: Type Codes
     *
     *  We store a 32-bit integer at the start of the buffer indicating
     *  which geometric type is being serialized.
     **********************************************************************/
    enum : std::uint32_t {
        TYPE_VEC  = 1,   // For std::vector< Vec >
        TYPE_AABB = 2,   // For std::vector< AABB<Vec> >
        TYPE_RAY  = 3,   // For std::vector< Ray<Vec> >
        TYPE_SEGM = 4,   // For std::vector< Segm<Vec> >
        TYPE_TRI  = 5    // For std::vector< Tri<Vec> >
    };

}
#endif //TYPECODES_H
