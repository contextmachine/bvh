//
// Created by Andrew Astakhov on 23.01.25.
//

#ifndef HELPERS_H
#define HELPERS_H
#include <vector>
#include "vec.h"
namespace bvh {
    inline void generateRayGrid(const vec3 &origin,
        const vec3 &direction,std::vector<Ray<vec3>> &grid,

        const double u=7.0,
        const double v=7.0,
        const int u_count=20,
        const int v_count=20){
        vec3 g_origin;
        grid.resize(u_count*v_count);
        for(int i=0; i<u_count; i++){
            for(int j=0; j<v_count; j++){
                grid[i*v_count+j] = {origin, g_origin+vec3(u*i, v*j, 0)-origin};
            }
        }
    }

}
#endif //HELPERS_H
