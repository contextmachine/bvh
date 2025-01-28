import numpy as np
from bvh import TriangleSoup, BVHTree

np.set_printoptions(suppress=True)

##############################################################################
# Example placeholders for user-provided data:
##############################################################################

# Replace the "..." placeholders with your own data initialization
mesh = TriangleSoup(...)        # Your triangle soup mesh
seats_bvh = BVHTree(...)        # Your BVH tree for "seats"
directions = ...                # Array-like of ray directions
starts = ...                    # Array-like of ray starting points

# The maximum number of reflections or bounces to process.
# If this was defined somewhere else in your code, bring that definition here.
max_order = 2

##############################################################################
# Utility / Data Structures
##############################################################################

class SimulationResult:
    """
    Stores the results of the ray simulation:
      - hits:      For each BVH leaf index, the intersection data (if any).
      - rays:      For each BVH leaf index, the accumulated ray segments.
      - dists:     For each BVH leaf index, the accumulated path lengths.
    """
    def __init__(self, hits, rays_path, rays_paths_length):
        self.hits = hits
        self.rays = rays_path
        self.dists = rays_paths_length

##############################################################################
# Helper function
##############################################################################

def record_hits(
    ray_paths_dict,         # dict: leaf_index -> list of lists of segment endpoints
    hits_dict,              # dict: leaf_index -> list of intersection hits
    path_lengths_dict,      # dict: leaf_index -> list of total travel distances
    raycast_index_dict,     # dictionary from seats_bvh.raycast()
    raycast_values,         # values returned by seats_bvh.raycast()
    global_indices,         # current global indices of the rays
    all_rays_cache,         # rrrays in original code, storing [start, direction]
    dists,                  # array of accumulated distances
    partial_paths_storage,  # rrays in original code: a list of lists of segments
    reflection_rays2        # the reflection array from mesh.reflection()
):
    """
    This function stores information about rays that actually hit a seat
    (based on the raycast result). The logic here is taken directly from the
    original dict_to_exists_csr, but renamed and cleaned up.
    """
    # raycast_index_dict is a dict with keys = (leaf_id, local_ray_index),
    # and values = something like an index in raycast_values array, etc.
    for (leaf_id, local_ray_idx), raycast_val_idx in raycast_index_dict.items():
        # Ensure the leaf has an entry in the dictionaries
        if leaf_id not in ray_paths_dict:
            ray_paths_dict[leaf_id] = []
            hits_dict[leaf_id] = []
            path_lengths_dict[leaf_id] = []

        # The original code checks the dot product to decide storing the hit
        # reflection_rays2[local_ray_idx, 0, :] => new start
        # all_rays_cache[global_indices[local_ray_idx]][0] => old start
        # all_rays_cache[global_indices[local_ray_idx]][1] => old direction
        #
        # raycast_values[raycast_val_idx, -1] => the t-parameter or distance factor
        old_start = all_rays_cache[global_indices[local_ray_idx]][0]
        old_dir   = all_rays_cache[global_indices[local_ray_idx]][1]
        new_start = reflection_rays2[local_ray_idx, 0, :]

        # Dist factor (last element in 'raycast_values')
        t_dist = raycast_values[raycast_val_idx, -1]

        # Check the dot product condition (same as original)
        dot_val = np.dot(new_start - old_start, old_dir)
        if dot_val >= t_dist:
            # Build the full path history (a list-of-lists structure)
            # partial_paths_storage[global_indices[local_ray_idx]] => the path so far
            # plus the new segment from old_start => old_start + old_dir * t_dist
            existing_path_list = partial_paths_storage[global_indices[local_ray_idx]]
            final_point = (old_start + old_dir * t_dist).tolist()

            # "ll" is the augmented path including the new segment
            new_path_history = np.array(existing_path_list).tolist() + [
                [old_start.tolist(), final_point]
            ]

            ray_paths_dict[leaf_id].append(new_path_history)
            path_lengths_dict[leaf_id].append(dists[global_indices[local_ray_idx]] + t_dist)
            hits_dict[leaf_id].append(raycast_values[raycast_val_idx, :-1])


##############################################################################
# Main Simulation Logic
##############################################################################

# Initialize the input rays array: shape (N, 2, 3),
# where rays[i,0] = start of i-th ray, rays[i,1] = direction of i-th ray.
rays = np.stack([np.array(starts), np.array(directions)], axis=1)

# Ensure the mesh has its BVH built
if not mesh.has_bvh():
    mesh.build_bvh()

# We keep rrrays (copy of the original rays) and rrays (list of partial paths)
# exactly as in the original code, for correct indexing/tracking.
all_rays_cache = np.copy(rays)        # rrrays
partial_paths_storage = [[] for _ in rays]  # rrays

# Accumulated distances per ray
dists = np.zeros(rays.shape[0], dtype=float)

# Each ray initially indexed by 0..N-1
indices_initial = np.arange(len(rays), dtype=int)
indices = np.copy(indices_initial)

# Dictionaries that will map leaf_id -> data
rays_paths = {}         # leaf_id -> list of path-segment histories
rays_hits = {}          # leaf_id -> list of hit data
rays_paths_lengths = {} # leaf_id -> list of path distances

for step in range(max_order):
    # Update all_rays_cache for the current "active" indices
    all_rays_cache[indices] = rays

    # Prepare arrays for reflection
    reflection_counts = np.empty(rays.shape[0], dtype=bool)
    reflection_rays2 = np.empty(rays.shape)

    # Perform reflection on the mesh (mesh must provide reflection(...) API)
    # reflection_counts[i] indicates if reflection was valid (e.g., if it actually hits)
    mesh.reflection(rays, reflection_rays2, reflection_counts)

    # Perform a raycast on seats_bvh (returns a dict of hits and an array of values)
    raycast_index_dict, raycast_values = seats_bvh.raycast(rays)

    # Keep track of which rays had a valid reflection
    valid_reflection_mask = reflection_counts > 0

    # Record hits (the function name is changed, but logic is identical)
    record_hits(
        ray_paths_dict=rays_paths,
        hits_dict=rays_hits,
        path_lengths_dict=rays_paths_lengths,
        raycast_index_dict=raycast_index_dict,
        raycast_values=raycast_values,
        global_indices=indices,
        all_rays_cache=all_rays_cache,
        dists=dists,
        partial_paths_storage=partial_paths_storage,
        reflection_rays2=reflection_rays2
    )

    # Update distances for the rays that validly reflected
    dists[indices[valid_reflection_mask]] += np.linalg.norm(
        rays[valid_reflection_mask, 0, :] - reflection_rays2[valid_reflection_mask, 0, :],
        axis=-1
    )

    # If no rays remain active, break out
    if not np.any(valid_reflection_mask):
        break

    # For valid reflections, the new direction's "start" is reflection_rays2[...,0]
    # For invalid ones, we just shoot a big distance in the old direction
    rays[valid_reflection_mask, 1, :] = reflection_rays2[valid_reflection_mask, 0, :]
    rays[~valid_reflection_mask, 1, :] = (
        rays[~valid_reflection_mask, 0, :] +
        rays[~valid_reflection_mask, 1, :] * 20000
    )

    # Store the updated rays in partial_paths_storage
    rays_list = rays.tolist()
    for local_i, global_i in enumerate(indices):
        partial_paths_storage[global_i].append(rays_list[local_i])

    # Now focus on only the still-active rays
    # reflection_rays2[valid_reflection_mask] -> shape (M, 2, 3)
    active_rays2 = reflection_rays2[valid_reflection_mask]

    # Normalize the direction vectors of active rays
    norm_factor = np.linalg.norm(active_rays2[..., 1, :], axis=1).reshape((-1, 1))
    active_rays2[..., 1, :] /= norm_factor

    # Update indices and rays for the next iteration
    indices = indices[valid_reflection_mask]
    rays = np.ascontiguousarray(active_rays2)

##############################################################################
# Post-processing: gather final data into lists
##############################################################################

all_leaf_hits = []
all_leaf_paths = []
all_leaf_lengths = []

leaf_count = seats_bvh.leafs_count()
for leaf_id in range(leaf_count):
    if leaf_id in rays_paths:
        all_leaf_paths.append(rays_paths[leaf_id])
        all_leaf_hits.append(np.array(rays_hits[leaf_id]).tolist())
        all_leaf_lengths.append(np.array(rays_paths_lengths[leaf_id]).tolist())
    else:
        all_leaf_paths.append([])
        all_leaf_hits.append([])
        all_leaf_lengths.append([])

##############################################################################
# Final result object
##############################################################################

simulation_result = SimulationResult(
    hits=all_leaf_hits,
    rays_path=all_leaf_paths,
    rays_paths_length=all_leaf_lengths
)

# 'simulation_result' now holds the final data 
# (equivalent to 'c' in the original code)