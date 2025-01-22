import json
import time

import numpy as np

from bvh import TriangleSoup
from pathlib import Path
inputs_triangles=Path(__file__).parent/"test_inputs_triangles.txt"
inputs_start=Path(__file__).parent/"test_inputs_start.txt"
inputs_direction=Path(__file__).parent/"test_inputs_direction.txt"
outputs_counts=Path(__file__).parent/"test_outputs_counts.txt"
outputs_pts=Path(__file__).parent/"test_outputs_pts.txt"
def read_array_from_path(path, dtype=None):
    with path.open('r') as f:
        return np.array(json.loads(f.read()),dtype=dtype)

def create_rays(starts, directions):
    return np.stack([starts, directions],
                    axis=1)
def test_raycast_hits():

    mesh=TriangleSoup(read_array_from_path(inputs_triangles))
    mesh.build_bvh()
    start=np.array(read_array_from_path(inputs_start))
    direction=read_array_from_path(inputs_direction)
    result=mesh.raycast_hits(create_rays(start,direction))
    hits,counts=np.array(result[0]),np.array(result[1],dtype=np.uint)
    out_counts=read_array_from_path(outputs_counts, dtype=np.uint)

    assert np.all(np.isclose(counts,out_counts))

    assert np.all(np.isclose(np.array(hits), read_array_from_path( outputs_pts)))

def test_raycast_counts():

    mesh=TriangleSoup(read_array_from_path(inputs_triangles))
    mesh.build_bvh()

    start=np.array(read_array_from_path(inputs_start))
    direction=read_array_from_path(inputs_direction)
    rays=create_rays(start, direction)
    counts = np.empty((rays.shape[0],), dtype=np.uint)
    mesh.raycast_counts(rays,counts)




    assert np.all(np.isclose(counts, read_array_from_path(outputs_counts,dtype=np.uint)))




if __name__ == '__main__':
    time.sleep(1.9)
    test_raycast_hits()
    test_raycast_counts()