//
// Created by Andrew Astakhov on 15.01.25.
//

#ifndef RAY_H
#define RAY_H

namespace bvh {
  template<typename Vec>
    struct Ray {
      Vec start;
      Vec dir;
    };

}
#endif //RAY_H
