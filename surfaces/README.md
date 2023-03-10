# Surface utilities

## fsLR-32k
| **File**                                | **Type** | **Description**                                  |
|-----------------------------------------|----------|--------------------------------------------------|
| fsLR-32k.?.sphere.reg.surf.gii          | gifti surface  | sphere to register between fsnative and fsLR-32k |
| fsLR-32k.?.inflated.surf.gii            | gifti surface  | Inflated canonical surface for data projection   |
| fsLR-32k.?.surf.gii                     | gifti surface  | Canonical surface for data projection            |

## fsLR-5k
| **File**                               | **Type** | **Description**                                  |
|----------------------------------------|----------|--------------------------------------------------|
| fsLR-5k.?.sphere.reg.surf.gii          | gifti surface  | ???    to register between fsnative and fsLR-32k |
| fsLR-5k.?.sphere.reg.deformed.surf.gii | gifti surface  | ???                                              |
| fsLR-5k.?.inflated.surf.gii            | gifti surface  | Inflated canonical surface for data projection   |
| fsLR-5k.?.surf.gii                     | gifti surface  | Canonical surface for data projection            |
| fsLR-5k.?.mask.shape.gii               | gifti data     | Medial wall mask int{0,1}                        |

## fsaverage5
| **File**                         | **Type** | **Description** |
|----------------------------------|----------|-----------------------|
| fsaverage5.?.sphere.reg.surf.gii | gifti surface  | sphere to register between fsnative and fsaverage5         |
| fsaverage.?.midthickness_orig.32k_fs_LR.surf.gii   | gifti surface  | ???           |

## fsaverage5/surf
| **File**      | **Type**         | **Description**                                                  |
|---------------|------------------|------------------------------------------------------------------|
| lh.pial       | freesurface surf | GM/CSF boundary surface                                          |
| lh.white      | freesurface surf | WM/GM boundary surface                                           |
| lh.inflated   | freesurface surf | smoothed and inflated to approach sphere, using fixed contraints |
| ?h.sphere     | freesurface surf | surface inflated all the way to a sphere                         |
| ?h.sphere.reg | freesurface surf | Registration sphere between native and fsaverage5 surf           |
