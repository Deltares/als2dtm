using StatsBase
using Statistics

# get 2D xy array from pointcloud
"read xy from Cloud"
function xy_array(cloud::Cloud, index::AbstractVector)
    n = length(index)
    XY = zeros(2, n)
    for i = 1:n
        p = positions(cloud)[index[i]]
        XY[1, i] = p[1]
        XY[2, i] = p[2]
    end
    XY
end
xy_array(cloud::Cloud) = xy_array(cloud, 1:length(cloud))

# get arrays of x/y/z coordinates
getx(cloud::Cloud, index::Vector{U}) where U <: Integer = getindex.(positions(cloud)[index], 1)
getx(cloud::Cloud, index::Integer) = positions(cloud)[index][1]
getx(cloud::Cloud) = getindex.(positions(cloud), 1)
gety(cloud::Cloud, index::Vector{U}) where U <: Integer = getindex.(positions(cloud)[index], 2)
gety(cloud::Cloud, index::Integer) = positions(cloud)[index][2]
gety(cloud::Cloud) = getindex.(positions(cloud), 2)
getz(cloud::Cloud, index::Vector{U}) where U <: Integer = getindex.(positions(cloud)[index], 3)
getz(cloud::Cloud, index::Integer) = positions(cloud)[index][3]
getz(cloud::Cloud) = getindex.(positions(cloud), 3)

# set arrays of x/y/z coordinates
setx!(cloud::Cloud, index::Vector{U}, value::Vector{Float64}) where U <: Integer = setindex!.(positions(cloud)[index], value, 1)
sety!(cloud::Cloud, index::Vector{U}, value::Vector{Float64}) where U <: Integer = setindex!.(positions(cloud)[index], value, 2)
setz!(cloud::Cloud, index::Vector{U}, value::Vector{Float64}) where U <: Integer = setindex!.(positions(cloud)[index], value, 3)
setx!(cloud::Cloud, index::Integer, value::Float64) = setindex!(positions(cloud)[index], value, 1)
sety!(cloud::Cloud, index::Integer, value::Float64) = setindex!(positions(cloud)[index], value, 2)
setz!(cloud::Cloud, index::Integer, value::Float64) = setindex!(positions(cloud)[index], value, 3)

# normalize z from given z value
function normalizez!(cloud::Cloud, index::Vector{U}, value::Float64) where U <: Integer
    z = getz(cloud, index)
    setz!(cloud, index, z - value)
end

# statistical functions
reducer_minz(cloud::Cloud, index::Vector{U}) where U <: Integer = minimum(getz(cloud, index))
reducer_maxz(cloud::Cloud, index::Vector{U}) where U <: Integer = maximum(getz(cloud, index))
reducer_medz(cloud::Cloud, index::Vector{U}) where U <: Integer = median(getz(cloud, index))
reducer_meanz(cloud::Cloud, index::Vector{U}) where U <: Integer = mean(getz(cloud, index))
reducer_madz(cloud::Cloud, index::Vector{U}) where U <: Integer = mad(getz(cloud, index), normalize=false)
reducer_stdz(cloud::Cloud, index::Vector{U}) where U <: Integer = std(getz(cloud, index))
reducer_99pz(cloud::Cloud, index::Vector{U}) where U <: Integer = maximum(trim(getz(cloud, index), prop=0.01))

# reducer for vegetation histogram splits
function reducer_veg_below_2m(cloud::Cloud, index::Vector{U}) where U <: Integer
    heights = getz(cloud, index)
    below = heights .< 2.
    return sum(below) / length(heights) * 100.
end
function reducer_veg_below_5m(cloud::Cloud, index::Vector{U}) where U <: Integer
    heights = getz(cloud, index)
    below = heights .< 5.
    return sum(below) / length(heights) * 100.
end
function reducer_veg_above_25m(cloud::Cloud, index::Vector{U}) where U <: Integer
    heights = getz(cloud, index)
    above = heights .> 25.
    return sum(above) / length(heights) * 100.
end
function reducer_veg_above_20m(cloud::Cloud, index::Vector{U}) where U <: Integer
    heights = getz(cloud, index)
    above = heights .> 20.
    return sum(above) / length(heights) * 100.
end

function reducer_ranz(cloud::Cloud, index::Vector{U}) where U <: Integer
    z = getz(cloud, index)
    maximum(z) - minimum(z)
end

function reducer_pctz(cloud::Cloud, index::Vector{U}) where U <: Integer
    quantile(getz(cloud, index), [0.1, 0.25, 0.5, 0.75, 0.9])
end

# reducer function to count number of points per cell
reducer_count(cloud::XYZ.Cloud, index::Vector{U}) where U <: Integer = length(index)
reducer_dens = reducer_count  # backwards compatibility

# reducer functions for any point attribute
function reducer_min_attr(field::Symbol)
    reducer(cloud::XYZ.Cloud, index::Vector{U}) where U <: Integer = minimum(cloud[field][index])
end
function reducer_max_attr(field::Symbol)
    reducer(cloud::XYZ.Cloud, index::Vector{U}) where U <: Integer = maximum(cloud[field][index])
end
function reducer_med_attr(field::Symbol)
    reducer(cloud::XYZ.Cloud, index::Vector{U}) where U <: Integer = median(cloud[field][index])
end
function reducer_min_abs_attr(field::Symbol)
    reducer(cloud::XYZ.Cloud, index::Vector{U}) where U <: Integer = minimum(abs.(cloud[field][index]))
end

# reducer function to get any attribute of the lowest or highest point
function reducer_lowestpoint_attr(field::Symbol)
    reducer(cloud::XYZ.Cloud, index::Vector{U}) where U <: Integer = cloud[field][index[argmin(getz(cloud, index))]]
end
function reducer_highestpoint_attr(field::Symbol)
    reducer(cloud::XYZ.Cloud, index::Vector{U}) where U <: Integer = cloud[field][index[argmax(getz(cloud, index))]]
end

# reduce_index, returns local(!) index of minimum / maximum
reducer_minz_index(cloud::Cloud, index::Vector{U}) where U <: Integer = argmin(getz(cloud, index))
reducer_maxz_index(cloud::Cloud, index::Vector{U}) where U <: Integer = argmax(getz(cloud, index))


const class_description = Dict{Int,String}(
    0 => "Created, never classified",
    1 => "Unclassified",
    2 => "Ground",
    3 => "Low Vegetation",
    4 => "Medium Vegetation",
    5 => "High Vegetation",
    6 => "Building",
    7 => "Low Point (noise)",
    8 => "Model Key-point (mass point)",
    9 => "Water",
    12 => "Overlap Points",

    18 => "High noise" # LAS 1.4
)

"Describe a `Cloud` with a number of statistics"
function describe(io::IO, cloud::Cloud)
    classes = sort!(unique(cloud[:classification]))
    quantvals = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    print(io, "quantiles,")
    join(io, quantvals, ",")
    println(io, "")
    for class in classes
        idxs = cloud[:classification] .== class
        println(io, "    class, $class, $(class_description[class])")
        print(io, "        elevation,")
        join(io, quantile(getz(cloud)[idxs], quantvals), ",")
        print(io, "\n        intensity,")
        join(io, quantile(cloud[:intensity][idxs], quantvals), ",")
        print(io, "\n")
    end
end

function describe(filepath::String, cloud::Cloud)
    open(filepath, "w") do io
        describe(io, cloud)
    end
end

const colorscale = 0x0100 # = 256

"Paint a Cloud efficiently with a Raster index and RGB bands from an orthophoto.
Does not load the entire orthophoto into memory, but one block at a time."
function paint!(pc::Cloud, r::Raster,
        rband::Ptr{GDAL.GDALRasterBandH},
        gband::Ptr{GDAL.GDALRasterBandH},
        bband::Ptr{GDAL.GDALRasterBandH})

    # 257 would be a full mapping, but 256 is in the standard,
    # see https://groups.google.com/d/msg/lasroom/Th2_-ywc2q8/4JRajnxxO0wJ
    # const colorscale = 0x0100 # = 256
    # assume all bands have same block sizes and data types
    @assert GDAL.getrasterdatatype(rband) === GDAL.GDT_Byte
    ncolblock, nrowblock = Ref(Cint(-1)), Ref(Cint(-1))
    GDAL.getblocksize(rband, ncolblock, nrowblock) # size of single blocks
    ncolblock, nrowblock = ncolblock[], nrowblock[] # get values from Ref
    ncol = GDAL.getrasterbandxsize(rband)
    nrow = GDAL.getrasterbandysize(rband)
    # number of blocks in the x and y direction
    nblockx = cld(ncol, ncolblock)
    nblocky = cld(nrow, nrowblock)
    # number of rows and columns in the last block of every row or column
    ncollastblock = rem(ncol, ncolblock) == 0 ? ncolblock : rem(ncol, ncolblock)
    nrowlastblock = rem(nrow, nrowblock) == 0 ? nrowblock : rem(nrow, nrowblock)
    # asserted it is GDAL.GDT_Byte
    # switched dimensions instead of transposing
    rblock = zeros(UInt8, ncolblock, nrowblock)
    gblock = zeros(UInt8, ncolblock, nrowblock)
    bblock = zeros(UInt8, ncolblock, nrowblock)
    # iterate over blocks
    for blockrow = 1:nblocky
        for blockcol = 1:nblockx
            nxblockoff = blockcol - 1
            nyblockoff = blockrow - 1
            # possible optimization: only read blocks that contain points?
            GDAL.readblock(rband, nxblockoff, nyblockoff, rblock)
            GDAL.readblock(gband, nxblockoff, nyblockoff, gblock)
            GDAL.readblock(bband, nxblockoff, nyblockoff, bblock)
            nrowthisblock = blockrow == nblocky ? nrowlastblock : nrowblock
            ncolthisblock = blockcol == nblockx ? ncollastblock : ncolblock
            # calculate the global indices for Cloud Raster index
            startrow = (blockrow - 1) * nrowblock
            startcol = (blockcol - 1) * ncolblock
            # iterate over block cells
            for i in 1:nrowthisblock
                for j in 1:ncolthisblock
                    row = startrow + i
                    col = startcol + j
                    idxs = r[row,col]
                    isempty(idxs) && continue # big speedup if many pixels don't have points
                    color = RGB{N0f16}(
                        reinterpret(N0f16, rblock[j,i] * colorscale),
                        reinterpret(N0f16, gblock[j,i] * colorscale),
                        reinterpret(N0f16, bblock[j,i] * colorscale))
                    pc[:color][idxs] = color
                end
            end
        end
    end
    nothing
end

function paint!(pc::Cloud, orthofile::AbstractString)
    # GDAL setup
    GDAL.allregister()
    ds = GDAL.open(orthofile, GDAL.GA_ReadOnly)
    geotransform = zeros(6)
    GDAL.getgeotransform(ds, geotransform)
    epsg = wkt2epsg(GDAL.getprojectionref(ds))
    rband = GDAL.getrasterband(ds, 1)
    gband = GDAL.getrasterband(ds, 2)
    bband = GDAL.getrasterband(ds, 3)
    ncol = Int(GDAL.getrasterbandxsize(rband))
    nrow = Int(GDAL.getrasterbandysize(rband))

    # create raster index that is guaranteed to align with the orthophoto pixels
    r = define_raster(pc, nrow, ncol, geotransform, epsg=epsg, pointfilter=nothing)
    paint!(pc, r, rband, gband, bband)
    GDAL.close(ds)
    nothing
end

function create_kdtree(cloud::Cloud;
    dim=2,
    metric=Euclidean(),
    index=1:length(cloud))

    # create tree
    if dim == 2
        tree = KDTree(xy_array(cloud, index), metric)
    elseif dim == 3
        tree = KDTree(positions(cloud)[index], metric)
    end

    tree
end

function point_density!(cloud::Cloud, spatialindex::NearestNeighbors.NNTree;
                        range=(1.0 / (4 / 3.0 * pi ))^(1.0 / 3),
                        pointfilter=nothing) # vol = 1m3
    @assert (spatialindex.metric == Euclidean()) && (length(spatialindex.data[1]) == 3)
    n = length(cloud)
    cloud[:point_density] = zeros(Float32, n)
    vol = 4 / 3.0 * pi * range^3

    # loop through points and calculate point density
    @showprogress 5 "looping through las points..." for i in 1:n
        (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
        # find number of points in range
        neighbors = inrange(spatialindex, positions(cloud)[i], range)
        # set point source attribute to visualize no. of neighbors
        cloud[:point_density][i] = length(neighbors) / vol
    end

    nothing
end
point_density!(cloud::Cloud; range=(1.0 / (4 / 3.0 * pi ))^(1.0 / 3)) = point_density!(cloud, create_kdtree(cloud; dim=3); range=range)

unique_id!(cloud::Cloud) = cloud[:id] = 1:length(cloud);

"""copy classification based on unique id,
assumes the field 'id' in cloud copy_from has the index of the same point in cloud copy_to"""
function copy_attributes!(cloud_copy_to::Cloud, cloud_copy_from::Cloud; field=:classification)
    @assert haskey(cloud_copy_from.attributes, :id) "both point clouds should have an 'id' field"
    @assert haskey(cloud_copy_from.attributes, field) "the 'copy_from' point cloud does not have a field $(field)"

    n = length(cloud_copy_to)
    if !haskey(cloud_copy_to.attributes, field)
        cloud_copy_to[field] = zeros(eltype(cloud_copy_from[field]), n)
    end

    # copy attribute valeu back to original cloud
    for (idx, v) in zip(cloud_copy_from[:id], cloud_copy_from[field])
        cloud_copy_to[field][idx] = v
    end
    nothing
    end

"""project the values of grid to point in pointcloud based on a raster index"""
function grid2cloud_attribute!(cloud::Cloud, grid::Array{T,2}, r::Raster, name::Symbol) where T <: Real
    @assert size(grid) == (r.nrow, r.ncol) "raster definition 'r' should have same size as 'grid'"
    cloud[name] = zeros(T, length(cloud))
    for icell in 1:length(r)
        v = grid[icell]
        for i in r[icell]
            cloud[name][i] = v
        end
    end
    nothing
end
