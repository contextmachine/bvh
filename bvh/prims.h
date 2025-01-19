//
// Created by Andrew Astakhov on 15.01.25.
//

#ifndef PRIMS_H
#define PRIMS_H


namespace bvh {
    template<typename Vec>
    /**
     * @struct Ray
     * @brief Represents a ray with a starting point and a direction.
     *
     * The Ray structure is primarily used in computational geometry, graphics, and physics simulations
     * to represent a line originating from a point and extending infinitely in a specific direction.
     *
     * @var start
     * The starting point of the ray, represented as a vector.
     *
     * @var direction
     * The direction of the ray, represented as a vector (unit or not unit).
     */
    struct Ray {
        Vec start;
        Vec direction;
    };

    template<typename Vec>
    /**
     * @struct Segm
     * @brief Represents a line segment defined by a start and an end point.
     *
     * The Segm structure is used in computational geometry to model a finite straight line
     * connecting two points in space.
     *
     * @var start
     * The starting point of the segment, represented as a vector.
     *
     * @var end
     * The ending point of the segment, represented as a vector.
     */
    struct Segm {
        Vec start;
        Vec end;
    };

    template<typename Vec>
    /**
     * @struct Tri
     * @brief Represents a triangle in a 3D space defined by its three vertices.
     *
     * The Tri structure is used in geometry processing, 3D rendering, and computational geometry
     * to represent a planar triangle formed by three points in space.
     *
     * @var a
     * The first vertex of the triangle, represented as a vector.
     *
     * @var b
     * The second vertex of the triangle, represented as a vector.
     *
     * @var c
     * The third vertex of the triangle, represented as a vector.
     */
    struct Tri {
        Vec a;
        Vec b;
        Vec c;
    };
}

#endif //PRIMS_H
