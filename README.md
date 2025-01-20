# BVH (Bounding Volume Hierarchy)

[![poetry-build](https://github.com/contextmachine/bvh/actions/workflows/poetry-build.yml/badge.svg)](https://github.com/contextmachine/bvh/actions/workflows/poetry-build.yml)

A high-performance C++/Cython implementation of Bounding Volume Hierarchy (BVH) for efficient ray-triangle intersection testing and spatial queries.

## Features

- Fast BVH construction for triangle meshes
- Efficient ray casting functionality
- Point-in-mesh testing
- Optimized C++ core with Python bindings
- Support for multiple Python versions (3.9+)

## Installation

The package can be installed using pip:

```bash
pip install git+https://github.com/contextmachine/bvh
```

For other platforms, you'll need to build from source, which requires:
- A C++ compiler
- Cython
- NumPy
- Poetry (for development)

## Usage

### Basic Example

```python
import numpy as np
from bvh import TriangleSoup

# Create a triangle mesh (n_triangles x 3 vertices x 3 coordinates)
triangles = np.array([
    [[0, 0, 0], [1, 0, 0], [0, 1, 0]],  # Triangle 1
    [[1, 0, 0], [1, 1, 0], [0, 1, 0]],  # Triangle 2
], dtype=np.float64)

# Initialize the BVH structure
soup = TriangleSoup(triangles)

# Build the BVH (required before any queries)
soup.build_bvh()

# Define rays for intersection testing
# Each ray is defined by origin and direction
rays = np.array([
    [[0, 0, 1], [1, 1, -1]],  # Ray 1: shooting down from above
], dtype=np.float64)

# Perform ray casting
hits, counts = soup.raycast_hits(rays)
```

### Point Inside Test

```python
# Test if points are inside the mesh
points = np.array([
    [0.5, 0.5, 0],  # Point to test
], dtype=np.float64)

inside = soup.points_inside(points)
```

## API Reference

### TriangleSoup

The main class for BVH operations.

#### Methods:

- `__init__(triangles)`: Initialize with a numpy array of triangles
- `build_bvh()`: Build the BVH structure
- `has_bvh()`: Check if BVH is built
- `raycast_counts(rays)`: Get intersection counts for rays
- `raycast_hits(rays)`: Get intersection points and counts
- `points_inside(points)`: Test if points are inside the mesh

## Development

The project uses Poetry for dependency management and building. To set up a development environment:

```bash
# Clone the repository
git clone https://github.com/contextmachine/bvh.git
cd bvh

# Install dependencies
poetry install

# Build the project
poetry build
```

## Testing

Tests are written using C++ and the test suite can be built using Meson:

```bash
meson setup buildDir
cd buildDir
meson test
```

## License

Licensed under the Apache License, Version 2.0 - see [LICENSE](LICENSE) for details.


## Technical Analysis

### Implementation Overview

This BVH implementation represents a practical balance between simplicity and functionality, with several notable characteristics:

#### Key Strengths

1. **Clean and Efficient Structure**
   - **Readable**:The code is easy to follow. The data structures are clearly defined, and the recursive building function follows the familiar pattern of "if leaf, store object; else split and recurse."
   - Array-based node layout using a single vector, improving memory locality
   - Simple but effective centroid-based splitting strategy
   - Memory-friendly implementation (2N-1 nodes for N objects)

2. **Practical Design Choices**
   - Straightforward top-down recursive build process
   - Splitting based on largest centroid extent - simpler than SAH but reasonably effective
   - Clean separation between splitting logic and tree construction
   - Focused scope: core BVH functionality without unnecessary complications

#### Future works
   - Could use median-of-medians approach instead of full sorting
   - Room for handling degenerate cases (overlapping centroids)
   - Could benefit from SAH splitting for more optimal trees
   - Opportunity for parallel build process in large scenes

### Best Use Cases
- Small to medium sized geometric scenes. 
    *The current implementation can compute up to 2 million rays per second on a single core (Apple M1-Pro) (from python code). This may be slightly faster from c++ code depending on the application, but not by orders of magnitude.*
- If you are looking for an implementation that you can easily adapt to your needs
- Projects requiring a clean, maintainable BVH implementation
- Educational contexts or reference implementations

