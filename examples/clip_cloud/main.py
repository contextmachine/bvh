from __future__ import annotations

import sys
from typing import TypedDict, Optional, Dict, NamedTuple

import numpy as np
from numpy.typing import NDArray
import pye57
from bvh import TriangleSoup
#/Volumes/CXM_GREY_01/yndx_les_points/ynx_les_e57/Yndx_00-31.e57
#pye57.E57('')
class E57ScanData(TypedDict):
    cartesianX: NDArray[float]
    cartesianY: NDArray[float]
    cartesianZ: NDArray[float]
    colorRed: Optional[NDArray[np.uint8]]
    colorGreen: Optional[NDArray[np.uint8]]
    colorBlue: Optional[NDArray[np.uint8]]
    intensity: Optional[NDArray[float]]




class E57Scan(NamedTuple):
    num: int
    header:pye57.ScanHeader
    scan:E57ScanData



def clip(pts:NDArray[float],mesh:TriangleSoup)->NDArray[int]:
    if not mesh.has_bvh():
        mesh.build_bvh()
    global_points= pts
    is_pt_in_bbox = np.zeros((global_points.shape[0],), bool)
    mesh.in_root_bbox(global_points, is_pt_in_bbox)

    pts_in_bbox = global_points[is_pt_in_bbox]
    if len(pts_in_bbox)==0:
        return np.array([],dtype=int)

    is_pt_in_mesh = np.zeros((pts_in_bbox.shape[0],), dtype=bool)
    mesh.points_inside(pts_in_bbox, is_pt_in_mesh)
    return np.arange(pts.shape[0],dtype=int)[is_pt_in_bbox][is_pt_in_mesh]

def clip_e57(mesh: TriangleSoup, source: pye57.E57 | str, target: pye57.E57 | str):
        if isinstance(source, str):
            source = pye57.E57(source)
        if isinstance(target, str):
            target = pye57.E57(target, 'w')



        cnt=source.scan_count
        for i in range(source.scan_count):
                print(f'{i}/{cnt}',end='\r',flush=True)
                header, scan = source.get_header(i), source.read_scan_raw(i)
                points = np.column_stack([scan['cartesianX'], scan['cartesianY'], scan['cartesianZ']])
                pts=source.to_global(points, header.rotation, header.translation)

                indices: NDArray[int] = clip(pts, mesh)
                if len(indices)==0:
                    continue
                new_scan = {k: v[indices] for k, v in scan.items()}

                target.write_scan_raw(new_scan, translation=header.translation, rotation=header.rotation)


        return True

if __name__ == '__main__':
    from pathlib import Path

    with open( Path(__file__).parent/'clipmesh.json') as f:
        tris=np.array(eval(f.read()))
    from bvh import TriangleSoup
    mesh=TriangleSoup(tris)
    mesh.build_bvh()
    e57 = pye57.E57('Yndx_031.e57')
    clip_e57(mesh, e57, 'Yndx_00-31-out.e57')
