"""
    trimesh_collapse(msh::DMesh{D,T,Simplex{2},3}; tol=1e-3) -> DMesh
    trimesh_collapse(msh::DMesh{D,T,Simplex{2},3}, nb; tol=1e-3) -> DMesh

Collapse short edges in a triangular mesh, returning a new mesh with fewer vertices and triangles.

Any edge with Euclidean length below `tol` is collapsed: one endpoint is removed and all
references to it are redirected to the surviving endpoint. When one endpoint is a boundary
vertex and the other is not, the boundary vertex is always kept. Otherwise the lower-index
vertex is kept.

For interior edges, both adjacent triangles are removed. For boundary edges, only the one
triangle containing the edge is removed. Neighboring triangles that become degenerate
(two or more vertices mapped to the same point) are also filtered out.

A pre-computed neighbor matrix `nb` from `element_face_neighbors` may be supplied to
avoid recomputing it.

# Arguments
- `msh`: Input triangular mesh (`Simplex{2}` elements, any spatial dimension)
- `nb`: Optional pre-computed element face neighbor matrix
- `tol`: Edge length threshold below which edges are collapsed (default: `1e-3`)

# Returns
- A new `DMesh` with short edges removed.

# Example
```julia
msh2 = trimesh_collapse(msh, tol=0.01)
```
"""
function trimesh_collapse(msh::DMesh{D,T,Simplex{2},3}; tol=1e-3) where {D,T}
    nb = element_face_neighbors(msh)
    trimesh_collapse(msh, nb; tol)
end

function trimesh_collapse(msh::DMesh{D,T,Simplex{2},3}, nb; tol=1e-3) where {D,T}
    np = length(msh.p)
    nt = length(msh.t)
    fmap = facemap(Simplex{2}())

    # --- Phase 1: Identify boundary vertices ---
    is_boundary = falses(np)
    for iel in 1:nt
        for jf in 1:3
            if nb[jf, iel][1] == 0
                edge = msh.t[iel][fmap[jf]]
                is_boundary[edge[1]] = true
                is_boundary[edge[2]] = true
            end
        end
    end

    # --- Phase 2: Mark collapses ---
    remap = collect(1:np)           # remap[i] = surviving vertex for i
    dead_tris = falses(nt)

    # Helper: follow remap chain to find current live representative
    function resolve(v)
        while remap[v] != v
            v = remap[v]
        end
        return v
    end

    for iel in 1:nt
        for jf in 1:3
            t2, _ = nb[jf, iel]

            # Skip interior edges where we'd process the same edge from t2 side
            t2 > 0 && t2 < iel && continue

            edge = msh.t[iel][fmap[jf]]
            va = resolve(edge[1])
            vb = resolve(edge[2])

            # Already collapsed to same point
            va == vb && continue

            # Check edge length
            norm(msh.p[va] - msh.p[vb]) >= tol && continue

            # Decide which vertex to keep: prefer boundary vertex
            v_keep, v_del = if is_boundary[va] && !is_boundary[vb]
                va, vb
            elseif is_boundary[vb] && !is_boundary[va]
                vb, va
            else
                min(va, vb), max(va, vb)
            end

            # Apply collapse
            remap[v_del] = v_keep

            # Mark triangles as dead
            dead_tris[iel] = true
            if t2 > 0
                dead_tris[t2] = true
            end
        end
    end

    # --- Phase 3: Path-compress remap chains (A→B→C becomes A→C) ---
    for i in 1:np
        remap[i] = resolve(i)
    end

    # --- Phase 4: Build compacted vertex list ---
    # Surviving vertices: those that map to themselves
    keep_verts = findall(i -> remap[i] == i, 1:np)

    # new_idx[old_index] = new index in compacted array (only valid for kept verts)
    new_idx = zeros(Int, np)
    for (k, v) in enumerate(keep_verts)
        new_idx[v] = k
    end

    # Composed map: old vertex i → new index
    final_idx(i) = new_idx[remap[i]]

    # --- Phase 5: Build compacted element list ---
    new_t = SVector{3, Int}[]
    sizehint!(new_t, nt - count(dead_tris))

    for iel in 1:nt
        dead_tris[iel] && continue

        tri = msh.t[iel]
        a, b, c = final_idx(tri[1]), final_idx(tri[2]), final_idx(tri[3])

        # Drop degenerate elements (two or more vertices collapsed to same point)
        (a == b || b == c || a == c) && continue

        push!(new_t, SVector(a, b, c))
    end

    return DMesh(msh.p[keep_verts], new_t)
end

"""
    trimesh_flip!(msh::DMesh{D,T,Simplex{2},3})
    trimesh_flip!(msh::DMesh{D,T,Simplex{2},3}, nb)

Improve triangular mesh quality by flipping edges between adjacent triangles.

Modifies the mesh in-place by examining each edge shared between two adjacent triangles
and flipping it across the opposite diagonal if the quality of both triangles improves
significantly. This operation is commonly used in mesh optimization algorithms.

The quality of a flip is evaluated using the mean ratio quality metric. An edge is flipped
if the minimum quality of the two new triangles exceeds the minimum quality of the original
pair by at least 0.025. For 3D meshes, the flip is additionally checked to ensure that
the resulting triangles have consistent normal orientations.

A pre-computed neighbor matrix `nb` from `element_face_neighbors` may be supplied to
avoid recomputing it.

# Arguments
- `msh`: Input triangular mesh (`Simplex{2}` elements, 2D or 3D)
- `nb`: Optional pre-computed element face neighbor matrix

# Example
```julia
trimesh_flip!(msh)  # Improve mesh quality in-place
```
"""
function trimesh_flip!(msh::DMesh{D,T,Simplex{2},3}) where {D,T}
    nb = element_face_neighbors(msh)
    trimesh_flip!(msh, nb)
end

function trimesh_flip!(msh::DMesh{D,T,Simplex{2},3}, nb) where {D,T}
    function trinormal3(pts)
        p1, p2, p3 = pts
        c = cross(p2 - p1, p3 - p1)
        return c / norm(c)
    end
    
    triqual3(pts) = element_quality_mean_ratio(Simplex{2}(), pts)

    nt = length(msh.t)
    for t1 in 1:nt
        for n1 in 1:3
            t2, n2 = nb[n1, t1]
            
            if t2 > 0
                old_t1 = msh.t[t1]
                old_t2 = msh.t[t2]
                
                q1 = triqual3(msh.p[old_t1])
                q2 = triqual3(msh.p[old_t2])
                minqold = min(q1, q2)
                
                if minqold < 0.9
                    tix11 = mod(n1, 3) + 1
                    tix12 = mod(n1 + 1, 3) + 1
                    tix21 = mod(n2, 3) + 1
                    tix22 = mod(n2 + 1, 3) + 1

                    # Create new triangles with swapped edges (no allocations with ntuple)
                    newt1 = SVector(ntuple(i -> i == tix12 ? old_t2[n2] : old_t1[i], 3))
                    newt2 = SVector(ntuple(i -> i == tix22 ? old_t1[n1] : old_t2[i], 3))

                    q3 = triqual3(msh.p[newt1])
                    q4 = triqual3(msh.p[newt2])
                    minqnew = min(q3, q4)

                    if minqnew > minqold + 0.025
                        flip = if D == 2
                            true  # In 2D, quality function handles orientation
                        else  # D == 3
                            normal1 = trinormal3(msh.p[old_t1])
                            normal2 = trinormal3(msh.p[old_t2])
                            normal3 = trinormal3(msh.p[newt1])
                            normal4 = trinormal3(msh.p[newt2])
                            if minqold < 0.01 # Desperate
                                (dot(normal3, normal4) > 0)
                            else
                                (dot(normal1, normal2) > 0) && (dot(normal3, normal4) > 0)
                            end
                        end
                        
                        if flip
                            # Update triangles in mesh
                            msh.t[t1] = newt1
                            msh.t[t2] = newt2

                            # Update neighbor connectivity
                            nbt, nbn = nb[tix21, t2]
                            nb[n1, t1] = (nbt, nbn)
                            if nbt > 0
                                nb[nbn, nbt] = (t1, n1)
                            end

                            nbt, nbn = nb[tix11, t1]
                            nb[n2, t2] = (nbt, nbn)
                            if nbt > 0
                                nb[nbn, nbt] = (t2, n2)
                            end

                            nb[tix11, t1] = (t2, tix21)
                            nb[tix21, t2] = (t1, tix11)
                        end
                    end
                end
            end
        end
    end
end
