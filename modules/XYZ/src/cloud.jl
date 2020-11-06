using Printf

"bounding box representing the smallest straightup rectangle containing all data "
struct BoundingBox
    xmin::Float64
    ymin::Float64
    zmin::Float64
    xmax::Float64
    ymax::Float64
    zmax::Float64
end

# construct BoundingBox with keyword argumets to prevent mixups
function BoundingBox(;xmin=0.0::Float64, ymin=0.0::Float64, zmin=0.0::Float64,
    xmax=0.0::Float64, ymax=0.0::Float64, zmax=0.0::Float64)
    BoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)
end

# #WARNING: Base.@sprintf is deprecated: it has been moved to the standard library package `Printf`.
# #Add `using Printf` to your imports.
# #https://docs.julialang.org/en/v1/stdlib/Printf/

# function Base.show(io::IO, bbox::BoundingBox)
#     print(io, @printf "xmin=%.1f, " bbox.xmin )
#     print(io, @printf "xmax=%.1f, " bbox.xmax )
#     print(io, @printf "ymin=%.1f, " bbox.ymin )
#     print(io, @printf "ymax=%.1f, " bbox.ymax )
#     print(io, @printf "zmin=%.2f, " bbox.zmin )
#     print(io, @printf "zmax=%.2f" bbox.zmax )
# end

"""Cloud is a point cloud container type similar to Clouds' Cloud, without
a spatial index, but with a simple bounding box that is kept up to date."""
mutable struct Cloud
    positions::Vector{MVector{3,Float64}}
    attributes::Dict{Symbol,Vector}
    boundingbox::BoundingBox
end

Cloud(pos, attrs) = Cloud(pos, attrs, calc_bbox(pos))

function Base.show(io::IO, cloud::Cloud)
    print(io, "PointCloud (N=$(length(cloud)), bbox=[$(cloud.boundingbox)], attributes=$(keys(cloud.attributes)))")
end


"Get the vector of point cloud positions"
positions(cloud::Cloud) = getfield(cloud, :positions)

"Get the Dict of point cloud attributes"
attributes(cloud::Cloud) = getfield(cloud, :attributes)

"Get the bounding box of the point cloud"
boundingbox(cloud::Cloud) = getfield(cloud, :boundingbox)

"Calculate bounding box of the positions"
function calc_bbox(P::Vector{SVector{3,Float64}})
    xmin, xmax, ymin, ymax, zmin, zmax = Inf, -Inf, Inf, -Inf, Inf, -Inf
    @inbounds for p in P
        xmin = min(xmin, p[1])
        xmax = max(xmax, p[1])
        ymin = min(ymin, p[2])
        ymax = max(ymax, p[2])
        zmin = min(zmin, p[3])
        zmax = max(zmax, p[3])
    end
    BoundingBox(xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax)
end

"Calculate bounding box of the positions"
function calc_bbox(P::Vector{MVector{3,Float64}})
    xmin, xmax, ymin, ymax, zmin, zmax = Inf, -Inf, Inf, -Inf, Inf, -Inf
    @inbounds for p in P
        xmin = min(xmin, p[1])
        xmax = max(xmax, p[1])
        ymin = min(ymin, p[2])
        ymax = max(ymax, p[2])
        zmin = min(zmin, p[3])
        zmax = max(zmax, p[3])
    end
    BoundingBox(xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax)
end

"Calculate bounding box of the point cloud"
function calc_bbox(cloud::Cloud)
    P = positions(cloud)
    calc_bbox(P)
end

"Update the bounding box of the point cloud"
function update_bbox!(cloud::Cloud)
    cloud.boundingbox = calc_bbox(cloud)
    nothing
end

"Subset a point cloud with specified indices"
function Base.getindex(cloud::Cloud, I)
    pos = positions(cloud)[I]
    attrs = Dict{Symbol,Vector}()
    for (k, v) in attributes(cloud)
        attrs[k] = v[I]
    end
    Cloud(pos, attrs)
        end

"Get a point cloud attribute by indexing with a `Symbol`"
Base.getindex(cloud::Cloud, attr::Symbol) = attributes(cloud)[attr]

"Set a point cloud attribute assigning to a `Symbol` index"
Base.setindex!(cloud::Cloud, A::AbstractVector, attr::Symbol) = attributes(cloud)[attr] = A

"Number of points in the point cloud"
Base.length(cloud::Cloud) = length(positions(cloud))
