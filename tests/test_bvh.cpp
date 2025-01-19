//
// Created by Andrew Astakhov on 08.12.24.
//
// Из корня репозитория
// clang++ -std=c++17 -o test/bvh_test test/bvh_test.cpp -march=armv8-a+simd -O3
// -I$PWD Или из папки ./test clang++ -std=c++17 -o bvh_test bvh_test.cpp
// -march=armv8-a+simd -O3 -I$PWD/..
#include <cassert>

#include "bvh/vec.h"
#include "bvh/bvh.h"
#include "bvh/aabb.h"
#include <chrono>
#include <cstring>
#include <iostream>

#include <vector>
#include "utils.h"
// Helper function to generate random float in range [min, max]

using namespace bvh;

int main(int argc, char* argv[]) {
  int num_objects = 10;  // Default number of objects is 100
  bool print_all = false; // Default is to not print all
  float space_size = 100.0f;

  // Parse command line arguments
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
      num_objects = std::atoi(argv[++i]); // Assign num_objects from the argument
    } else if (strcmp(argv[i], "-print-all") == 0) {
      print_all = true; // Set print_all flag if option is present
    }
  }

  try {
    // Define a few 3D objects with their AABBs

    // Define a target AABB for querying the BVH tree



    std::vector<AABB<double, 3>>objects;


    for (int i = 0; i < num_objects; ++i) {
      vec3 point1=random_vec3(0.,10.);
      vec3 point2=point1-random_vec3(1.,10.);
      vec3 min_point=point1.min(point2);
      vec3 max_point=point1.max(point2);
      objects.push_back(AABB<double, 3>(min_point, max_point));

    }
    if (print_all){
    std::cout<<"[";
    for (int i = 0; i < num_objects; ++i) {
      auto& obj=objects[i];
      std::cout<<obj;
      if (i<(num_objects-1)) {
        std::cout<<",";
      }
    }
    std::cout<<"]\n\n";
    }
    // Create vectors of pointers for each BVH construction



    const auto start = std::chrono::high_resolution_clock::now();
    auto three=BVH();
    three.build(objects);
    std::cout<<three.nodes[0].left<< "  "<<three.nodes[0].right<<std::endl;
    assert(three.nodes[0].left==1);
    assert(three.nodes[0].right==2);
    std::cout<<three<<std::endl;
    
    const auto end = std::chrono::high_resolution_clock::now();

    const auto duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                             start)
            .count();
    if (print_all) {
      std::cout << "\n\nBVH\n--------\n" << std::endl;

    }
    //std::cout  << tree << std::endl;
    if (print_all) {
      std::cout << "BVH build time:  \n    " << duration << " ns. \n    "<< duration*1e-9 << " secs."<< std::endl;
    }
    } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown error occurred" << std::endl;
    }
    return 0;
}