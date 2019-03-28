#=
functions to create index points from pointcloud along line defined by start & end point.

profiler generetes DataFrame with XYZ and attribute information from pointcloud,
where X and Y are defined along line

=#

struct Profile{U<:Integer,T<:AbstractFloat}
    pointindex::Vector{U}
    pointlx::Vector{T}
    pointly::Vector{T}
    maxdist::T
    pstart::Tuple{T, T}
    pend::Tuple{T, T}
    len::T
end

# length of profile line segment
linelength(p::Profile) = p.len
pointindex(p::Profile) = p.pointindex
pointlx(p::Profile) = p.pointlx
pointly(p::Profile) = p.pointly

# Base functions for Profile type
Base.length(p::Profile) = length(p.pointindex)

"index points along profileline"
function define_profile(cloud::Cloud, ps::Tuple{T,T}, pe::Tuple{T,T}, maxdist::T;
    pointfilter = nothing) where T <: Real

    # define line and bbox
    x_min, x_max = minmax(ps[1], pe[1])
    y_min, y_max = minmax(ps[2], pe[2])
    x_min -= maxdist
    y_min -= maxdist
    x_max += maxdist
    y_max += maxdist

    # define rotated grid orientation
    transrot = segment_transformation(ps, pe)
    len = hypot(ps[1] - pe[1], ps[2] - pe[2]) #length(ps, pe)

    # loop through Cloud and store points
    n = length(positions(cloud))
    idx = Int32[]  # index of points in Cloud
    lx = Float64[]  # distance along line
    ly = Float64[]  # distance from line
    for i in 1:n
        xi, yi, zi = positions(cloud)[i]
        outside_bbox(x_min, y_min, x_max, y_max, xi, yi) && continue
        (pointfilter != nothing) && (pointfilter(cloud, i) || continue)

        p = SVector(xi, yi)
        p′ = transrot(p)

        if (0.0 <= p′[1] < len) && (abs(p′[2]) <= maxdist)
            push!(idx, i)
            push!(lx, p′[1])
            push!(ly, p′[2])
        end
    end
    np = length(idx)
    # sort points based on distance along line lx
    if np >= 100_000 # RadixSort only faster with many points
        perm = sortperm(lx, alg=RadixSort)
    else
        perm = sortperm(lx)
    end

    Profile(idx[perm], lx[perm], ly[perm], maxdist, ps, pe, len)
end


"create new profile based on pointfilter applied to points in profile"
function filter_profile(cloud::Cloud, prof::Profile, pointfilter)
    n = length(prof)
    idx = Int32[]
    for i in 1:n
        idxi = pointindex(prof)[i]
        (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
        push!(idx, idxi)
    end

    Profile(
        pointindex(prof)[idx],
        pointlx(prof)[idx],
        pointly(prof)[idx],
        prof.maxwidth,
        prof.pstart,
        prof.pend,
        prof.len)
end

# utils
"""Set up a transformation that will map points to a rotated
and translated coordinate system with the origin on the segment
start, and positive x along the segment."""
function segment_transformation(ps, pe)
    θ = -atan2(pe[2]-ps[2], pe[1]-ps[1])
    rot =  SMatrix{2,2}([cos(θ) -sin(θ); sin(θ) cos(θ)])
    trans = SVector(-ps[1], -ps[2])
    transrot(x::SVector{2,Float64}) = *(rot, x+trans)
end


"check if a point falls outside a bounding box"
function outside_bbox(x_min, y_min, x_max, y_max, x, y)
    x < x_min || y < y_min || x > x_max || y > y_max
end
