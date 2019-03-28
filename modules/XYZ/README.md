# XYZ

Julia package for easily setting up pipelines for processing point clouds.
The goal is to have a simple API that is flexible enough for diverse processing, including processes
that take into account neighbouring points. Such a mockup API is presented below.

## High level API
```julia
using XYZ

pipeline(infile, process)
```

## Medium level API
```julia
using XYZ

pts = read("las.las")

pts::PointCloud

function predicate(p)
    true
end

function classify(p) # or colorize
    p′
end

filter!(pts, pointfilter)
classify!(pts, predicate, class)

laspts = (header, Vector{LasPoint}, spatialindex)
laspt = (header, Vector{LasPoint}[i], spatialindex)



custom_env_filter1!(pts)
custom_env_filter2!(pts)

write_raster(pts::Vector{LasPoint}) = write_raster(xyz)

write("out.las", pts::PointCloud, predicate, option)
write("out.tif", pts, predicate, option) # 2D grid
write("out.csv", pts, predicate, option) # 2D profile
write("out.csv", pts, predicate, option) # XYZ with extra columns


function PointCloud{T <: AbstractFloat}(laspoints::Vector{LasPoint}, header::LasHeader, )
    positions = StaticArray
    Dict{Symbol,Vector{Any}}()
    PointCloud(positions, BallTree(xy_array(header, lasp), Chebyshev()), attr, header)

end

cloud = PointCloud(positions)

cloud = PointCloud(positions, Chebyshev())
cloud = PointCloud(positions, Eucldean())

type PointCloud{Dim,T,SIndex}
    positions::Vector{SVector{Dim,T}}
    spatial_index::SIndex
    attributes::Dict{Symbol,Vector{Any}}
    header::LasHeader
end
```
