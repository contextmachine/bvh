//
// Created by Andrew Astakhov on 15.01.25.
//

#ifndef PRIMS_H
#define PRIMS_H


namespace bvh {
    template<typename Vec>
    struct Ray {
        Vec start;
        Vec direction;
    };

    template<typename Vec>
    struct Segm {
        Vec start;
        Vec end;
    };

    template<typename Vec>
    struct Tri {
        Vec a;
        Vec b;
        Vec c;
    };
}

#endif //PRIMS_H
