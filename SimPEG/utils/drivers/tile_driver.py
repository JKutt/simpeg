import numpy as np
from discretize.utils import mesh_builder_xyz, refine_tree_xyz, active_from_xyz
from discretize import TreeMesh
from SimPEG.maps import TileMap
from scipy.spatial import Delaunay, cKDTree


def create_tile_meshes(
    locations,
    topography,
    indices,
    base_mesh=None,
    max_distance=np.inf,
    core_cells=[10, 10, 10],
    locations_refinement=[5, 5, 5],
    topography_refinement=[0, 0, 2],
    padding_distance=[[0, 0], [0, 0], [0, 0]],
):

    assert isinstance(indices, list), "'indices' must be a list of integers"

    if len(padding_distance) == 1:
        padding_distance = [padding_distance*2]*3

    global_mesh = mesh_builder_xyz(
        locations, core_cells,
        padding_distance=padding_distance,
        mesh_type='TREE', base_mesh=base_mesh,
        depth_core=1000
    )
    local_meshes = []
    for ind in indices:
        local_mesh = mesh_builder_xyz(
            locations, core_cells,
            padding_distance=padding_distance,
            mesh_type='TREE', base_mesh=base_mesh,
            depth_core=1000
        )
        local_mesh = refine_tree_xyz(
            local_mesh, topography,
            method='surface', octree_levels=topography_refinement,
            finalize=False
        )
        local_mesh = refine_tree_xyz(
            local_mesh, locations[ind],
            method='surface', octree_levels=locations_refinement,
            max_distance=max_distance,
            finalize=True
        )
        global_mesh.insert_cells(
            local_mesh.gridCC,
            local_mesh.cell_levels_by_index(np.arange(local_mesh.nC)),
            finalize=False,
        )

        local_meshes.append(local_mesh)

    global_mesh.finalize()
    global_active = active_from_xyz(global_mesh, topography, method='linear')

    # Cycle back to all local meshes and create tile maps
    local_maps = []
    for mesh in local_meshes:
        local_maps.append(
            TileMap(global_mesh, global_active, mesh)
        )

    return (global_mesh, global_active), (local_meshes, local_maps)


def create_nested_mesh(
        locations,
        base_mesh,
        method="convex_hull",
        max_distance=100.,
        pad_distance=1000.,
        min_level=2,
        finalize=True
):
    nested_mesh = TreeMesh(
        [base_mesh.h[0], base_mesh.h[1], base_mesh.h[2]], x0=base_mesh.x0
    )

    min_level = (base_mesh.max_level - min_level)
    base_refinement = base_mesh.cell_levels_by_index(np.arange(base_mesh.nC))
    base_refinement[base_refinement > min_level] = min_level

    nested_mesh.insert_cells(
        base_mesh.gridCC,
        base_refinement,
        finalize=False,
    )

    tree = cKDTree(locations[:, :2])
    rad, _ = tree.query(base_mesh.gridCC[:, :2])

    indices = np.where(rad < pad_distance)[0]
    # indices = np.where(tri2D.find_simplex(base_mesh.gridCC[:, :2]) != -1)[0]
    levels = base_mesh.cell_levels_by_index(indices)
    levels[levels==base_mesh.max_level] = base_mesh.max_level-1
    nested_mesh.insert_cells(
        base_mesh.gridCC[indices, :],
        levels,
        finalize=False,
    )

    if method == "convex_hull":
        # Find cells inside the data extant
        tri2D = Delaunay(locations[:, :2])
        indices = tri2D.find_simplex(base_mesh.gridCC[:, :2]) != -1
    else:
        # tree = cKDTree(locations[:, :2])
        # rad, _ = tree.query(base_mesh.gridCC[:, :2])
        indices = rad < max_distance

    nested_mesh.insert_cells(
        base_mesh.gridCC[indices, :],
        base_mesh.cell_levels_by_index(np.where(indices)[0]),
        finalize=finalize,
    )

    return nested_mesh