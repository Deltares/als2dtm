using Nullables

#= 
Functions to create a raster from Clouds.

reg_grid and reg_grid_index create a link between raster cells and Points in the
Cloud.
rasterize summarizes the Points per cell to one statistical value with and creates
the actual rasters

Dirk Eilander, Deltares, 11-2016

CHANGELOG
Version 0.3f
-grid index is calculated based on coordinates instead of using KDTree (faster)
-eliminated use of las header. the bbox is derived from spatial index pointcloud
Version 0.2
-fixed points on boundaries (top en left are included, right and bottom excluded)
-created Raster type
-added pointfilter to define_raster & reg_grid_index
Version 0.1
-bugfix for orientation grid

TODO:
- =#

struct Raster{U <: Integer}
    pointindex::Vector{Vector{U}}
    nrow::Int
    ncol::Int
    bbox::BoundingBox
    cellsize::Float64
    epsg::Nullable{Int}
    npoints::Int # number of unique points in index
    mask::BitArray{2}  # mask == True -> missing
end

function Raster(
    pointindex::Vector{Vector{U}},
    nrow::Int,
    ncol::Int,
    bbox::BoundingBox,
    cellsize::Float64,
    epsg::Nullable{Int},
    npoints::Int
    ) where U <: Integer

    Raster(pointindex, nrow, ncol, bbox, cellsize, epsg, npoints, falses(nrow, ncol))
end

function Raster(
    pointindex::Vector{Vector{U}},
    nrow::Int,
    ncol::Int,
    bbox::BoundingBox,
    cellsize::Float64,
    epsg::Int,
    npoints::Int,
    mask::BitArray{2}
    ) where U <: Integer

    Raster(pointindex, nrow, ncol, bbox, cellsize, Nullable{Int}(epsg), npoints, mask)
end

function Base.show(io::IO, r::Raster)
    println(io, "RasterIndex (nrow=$(r.nrow), ncol=$(r.ncol), cellsize=$(r.cellsize), n unique points=$(r.npoints), epsg=$(r.epsg), bbox=[$(r.bbox)])")
end

# construct coordinates in vector or meshgrid from Raster
epsg(r::Raster) = r.epsg
npoints(r::Raster) = r.npoints # number of unique points in index

# length is number of cells
Base.length(r::Raster) = length(r.pointindex)
Base.size(r::Raster) = (r.nrow, r.ncol)
Base.eachindex(r::Raster) = eachindex(r.pointindex)

function ind2sub(a, i)
    i2s = CartesianIndices(a)
    i2s[i].I
end
function sub2ind(a, i...)
    s2i = LinearIndices(a)
    s2i[i...]
end

ind2sub(r::Raster, i::Integer) = ind2sub((r.nrow, r.ncol), i)
function sub2ind(r::Raster, row::Integer, col::Integer)
    @info "raster" (r.nrow, r.ncol)
    @info "req" (row, col)
    sub2ind((r.nrow, r.ncol), row, col)
end

"transform coordinates to row, col in raster"
function coords2sub(r::Raster, x::AbstractFloat, y::AbstractFloat)
    # cells include top en left boundaries
    col = Int(fld(x - r.bbox.xmin, r.cellsize) + 1)
    row = Int(fld(r.bbox.ymax - y, r.cellsize) + 1)
    # include points on lower and right bbox boundaries
    if x == r.bbox.xmax
      col = col - 1
    end
    if y == r.bbox.ymin
      row = row - 1
    end
    ((col <= 0) || (row <= 0) || (row > r.nrow) || (col > r.ncol)) && @error("coordinates outside raster domain")
    row, col
end

# based on index of raster pointindex
Base.getindex(r::Raster, i::Integer) = r.pointindex[i]
# based on index of raster, return row, col, pointindex
Base.getindex(r::Raster, row::Integer, col::Integer) = r[sub2ind(r, row, col)]
# based on coordinates raster, return pointindex
function Base.getindex(r::Raster, x::AbstractFloat, y::AbstractFloat)
    row, col = coords2sub(r, x, y)
    r[row, col]
end

"transform coordinates to row, col in raster"
function coords2ind(bbox::BoundingBox,
    nrow::Integer,
    ncol::Integer,
    cellsize::Float64,
    x::AbstractFloat,
    y::AbstractFloat)

    # (bbox.xmin < p.X <= bbox.xmax && bbox.ymin < p.Y <= bbox.ymax ) || continue # && u_min_z <= p.Z <= u_max_z) || continue
    # cells include bottom and left boundaries
    col = Int(fld(x - bbox.xmin, cellsize) + 1)
    row = Int(fld(bbox.ymax - y, cellsize) + 1)
    # include points on lower and right bbox boundaries
    if x == bbox.xmax
      col = col - 1
    end
    if y == bbox.ymin
      row = row - 1
    end
    isinside = !((col <= 0) || (row <= 0) || (row > nrow) || (col > ncol))
    if !isinside
        row, col = 1, 1
    end
    sub2ind((nrow, ncol), row, col), isinside
end

"transform coordinates to row, col in raster"
function coords2ind_overlap(bbox::BoundingBox,
    nrow::Integer,
    ncol::Integer,
    cellsize::Float64,
    x::AbstractFloat,
    y::AbstractFloat;
    overlap=100.0)

    # cells include bottom and left boundaries
    col = Int(fld(x - bbox.xmin, cellsize) + 1)
    row = Int(fld(bbox.ymax - y, cellsize) + 1)
    # include points on lower and right bbox boundaries
    # if x == bbox.xmax
    #   col = col - 1
    # end
    # if y == bbox.ymin
    #   row = row - 1
    # end
    isinside = !((col <= 0) || (row <= 0) || (row > nrow) || (col > ncol))
    # return a list of indices

    # this assumes that the overlap is less than half the cell size
    inleft = rem(x - bbox.xmin, cellsize) <= overlap
    inbottom = rem(y - bbox.ymin, cellsize) <= overlap
    inright = rem(x - bbox.xmin, cellsize) >= (cellsize - overlap)
    intop = rem(y - bbox.ymin, cellsize) >= (cellsize - overlap)

    # correct for boundary issues
    if !isinside
        inleft = false
        inbottom = false
        inright = false
        intop = false
    end
    if col == 1
        inleft = false
    end
    if row == nrow
        inbottom = false
    end
    if col == ncol
        inright = false
    end
    if row == 1
        intop = false
    end

    row, col, isinside, inleft, inbottom, inright, intop
end

"transform row, col in raster to center cell x, y coordinates"
function ind2coords(bbox::BoundingBox,
    nrow::Integer,
    ncol::Integer,
    cellsize::Float64,
    icell::Integer)

    row, col = ind2sub((nrow, ncol), icell)
    # cells include bottom and left boundaries
    x = bbox.xmin + (col - 0.5) * cellsize
    y = bbox.ymax - (row - 0.5) * cellsize

    x, y
end

function rastercellcoordinates(r::Raster)
    ncells = r.ncol * r.nrow
    cellcoords = Vector{SVector{2,Float64}}()
    sizehint!(cellcoords, ncells)
    for icell in 1:ncells
        coords = SVector(ind2coords(r.bbox, r.nrow, r.ncol, r.cellsize, icell))
        push!(cellcoords, coords)
    end
    cellcoords
end

"raster definition constructor from Cloud, cellsize and predefined tiling and epsg"
function define_raster(cloud::Cloud, orig_tile::BoundingBox, overlap::Float64, cellsize::Float64;
    snapgrid=cellsize,
    epsg=Nullable{Int}(),
    pointfilter=nothing)
    # determine inner col and rows and its outer col and rows
    ocol = round(Int, abs(orig_tile.xmax - orig_tile.xmin) / cellsize)
    orow = round(Int, abs(orig_tile.ymax - orig_tile.ymin) / cellsize)
    overlapcells = Int(cld(overlap, cellsize))
    ncol = round(Int, ocol + 2 * overlapcells)
    nrow = round(Int, orow + 2 * overlapcells)

    # create mask and mask out overlaps
    mask = trues(nrow, ncol)
    mask[(overlapcells + 1):(end - overlapcells), (overlapcells + 1):(end - overlapcells)] .= false

    # calculate bounding box including overlap
    bbox = BoundingBox(xmin=orig_tile.xmin - overlapcells * cellsize,
                       xmax=orig_tile.xmax + overlapcells * cellsize,
                       ymin=orig_tile.ymin - overlapcells * cellsize,
                       ymax=orig_tile.ymax + overlapcells * cellsize
                       )

    # create grid definition
    xy_idx, npoints = reg_grid_index(cloud, nrow, ncol, bbox, cellsize;
        pointfilter=pointfilter)
    Raster(xy_idx, nrow, ncol, bbox, cellsize, epsg, npoints, mask)
end

"raster definition constructor from Cloud, cellsize and epsg"
function define_raster(cloud::Cloud, cellsize::Float64;
    snapgrid=cellsize,
    epsg=Nullable{Int}(),
    pointfilter=nothing)

    # create grid definition
    nrow, ncol, bbox = reg_grid(cloud, cellsize; snapgrid=snapgrid)
    xy_idx, npoints = reg_grid_index(cloud, nrow, ncol, bbox, cellsize;
        pointfilter=pointfilter)
    Raster(xy_idx, nrow, ncol, bbox, cellsize, epsg, npoints)
end

"raster definition constructor from Cloud, and GDAL geotransform"
function define_raster(cloud::Cloud, nrow::Integer, ncol::Integer, geotransform::Vector{Float64};
    epsg=Nullable{Int}(),
    pointfilter=nothing)

    # create grid definition
    cellsize = Float64(geotransform[2])
    nrow, ncol, bbox = reg_grid(nrow, ncol, geotransform)
    xy_idx, npoints = reg_grid_index(cloud, nrow, ncol, bbox, cellsize;
        pointfilter=pointfilter)
    Raster(xy_idx, nrow, ncol, bbox, cellsize, epsg, npoints)
end

"raster definition constructor from Cloud, and GDAL geotransform"
function define_raster_overlap(cloud::Cloud, nrow::Integer, ncol::Integer, geotransform::Vector{Float64};
    epsg=Nullable{Int}(),
    pointfilter=nothing,
    overlap=100.0)

    # create grid definition
    cellsize = Float64(geotransform[2])
    nrow, ncol, bbox = reg_grid(nrow, ncol, geotransform)
    xy_idx, npoints = reg_grid_index_overlap(cloud, nrow, ncol, bbox, cellsize;
        pointfilter=pointfilter,
        overlap=overlap)
    Raster(xy_idx, nrow, ncol, bbox, cellsize, epsg, npoints)
end

"create vectors with bbox from Cloud spatial index based on cell size
snap to utm coordinates by default"
function reg_grid(cloud::Cloud, cellsize::Float64; snapgrid=cellsize)
    # snap bbox to cellsize coordinates (to the outside)
    # TODO this assumes the coordinates are in cm precision
    # do the fld with scaled integers to prevent floating point issues
    # like fld(549.0, 0.01) == 54899.0
    if snapgrid != 0
        snapgrid_int = round(Int32, 100 * snapgrid)
        x_min_int = round(Int32, 100 * boundingbox(cloud).xmin)
        y_min_int = round(Int32, 100 * boundingbox(cloud).ymin)
        x_max_int = round(Int32, 100 * boundingbox(cloud).xmax)
        y_max_int = round(Int32, 100 * boundingbox(cloud).ymax)
        x_min = fld(x_min_int, snapgrid_int) * snapgrid_int * 0.01
        y_min = fld(y_min_int, snapgrid_int) * snapgrid_int * 0.01
        x_max = cld(x_max_int, snapgrid_int) * snapgrid_int * 0.01
        y_max = cld(y_max_int, snapgrid_int) * snapgrid_int * 0.01
    end

    # determine size of grid
    dy = y_max - y_min
    dx = x_max - x_min
    nrow = round(Int, dy / cellsize)
    ncol = round(Int, dx / cellsize)
    bbox = BoundingBox(xmin=x_min, ymin=y_min, xmax=x_max, ymax=y_max, zmin=NaN, zmax=NaN)
    nrow, ncol, bbox
end

"create vectors with bbox based on a GDAL geotransform vector"
function reg_grid(nrow::Integer, ncol::Integer, geotransform::Vector{Float64})
    # assumes a square grid
    cellsize = Float64(geotransform[2])
    x_min = geotransform[1]
    y_max = geotransform[4]
    x_max = x_min + nrow * cellsize
    y_min = y_max - ncol * cellsize
    bbox = BoundingBox(xmin=x_min, ymin=y_min, xmax=x_max, ymax=y_max)

    # determine size of grid
    dy = y_max - y_min
    dx = x_max - x_min
    nrow = round(Int, dy / cellsize)
    ncol = round(Int, dx / cellsize)

    nrow, ncol, bbox
end

"Create reg_grid_index from Cloud and vectors of x & y coordinates"
function reg_grid_index(
    cloud::Cloud,
    nrow::Int,
    ncol::Int,
    bbox::BoundingBox,
    cellsize::Float64;
    pointfilter=nothing)

    n = nrow * ncol # number of cells
    # the idx array gives per LAS point the linear index into the grid it was assigned to
    idx = Array{Vector{Int32}}(undef, n) # some stay zero because of skipped points
    for icell in 1:n
        idx[icell] = Int32[]
    end

    npoints = 0
    @showprogress 2 "Indexing grid..."  for i in 1:length(cloud)
        (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
        # calc cellindex
        pos = positions(cloud)[i]
        xi, yi = pos[1], pos[2]
        icell, isinside = coords2ind(bbox, nrow, ncol, cellsize, xi, yi)
        isinside && push!(idx[icell], i)
        npoints += 1
    end

    idx, npoints
end

"Create reg_grid_index from Cloud and vectors of x & y coordinates"
function reg_grid_index_overlap(
    cloud::Cloud,
    nrow::Int,
    ncol::Int,
    bbox::BoundingBox,
    cellsize::Float64;
    pointfilter=nothing,
    overlap=100.0)

    n = nrow * ncol # number of cells
    # the idx array gives per LAS point the linear index into the grid it was assigned to
    idx = Array{Vector{Int32}}(n) # some stay zero because of skipped points
    for icell in 1:n
        idx[icell] = Int32[]
    end

    # with overlap the main difference is that each point may not be in 0 or 1 cells,
    # but 0 to n (normally n is max 4)
    npoints = 0
    for i in 1:length(cloud)
        (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
        # calc cellindex
        pos = positions(cloud)[i]
        xi, yi = pos[1], pos[2]
        row, col, isinside, inleft, inbottom, inright, intop = coords2ind_overlap(bbox, nrow, ncol, cellsize, xi, yi, overlap=overlap)
        icell = sub2ind((nrow, ncol), row, col)

        ileft = sub2ind((nrow, ncol), row, col - 1)
        ibottom = sub2ind((nrow, ncol), row + 1, col)
        iright = sub2ind((nrow, ncol), row, col + 1)
        itop = sub2ind((nrow, ncol), row - 1, col)

        ileftbottom = sub2ind((nrow, ncol), row + 1, col - 1)
        irightbottom = sub2ind((nrow, ncol), row + 1, col + 1)
        irighttop = sub2ind((nrow, ncol), row - 1, col + 1)
        ilefttop = sub2ind((nrow, ncol), row - 1, col - 1)

        # center
        isinside && push!(idx[icell], i)
        inleft && push!(idx[ileft], i)
        inbottom && push!(idx[ibottom], i)
        inright && push!(idx[iright], i)
        intop && push!(idx[itop], i)

        # don't only add overlap to the sides but also away from the corners
        inleft && inbottom && push!(idx[ileftbottom], i)
        inright && inbottom && push!(idx[irightbottom], i)
        inright && intop && push!(idx[irighttop], i)
        inleft && intop && push!(idx[ilefttop], i)

        # TODO not sure what this should be in case of overlap
        # and if it matters
        npoints += 1
    end

    idx, npoints
end
        
"""fill empty raster index cells with nearest neighbor point
this method creates a new raster index"""
function inpaint_missings_nn(r::Raster, cloud::Cloud;
    fill_max_dist=10)
    # index points that satisfy pointfilter
    index = Int32[]
    for i = 1:length(r)
        for idx in r[i]
            push!(index, idx)
        end
    end

    missings = falses((r.nrow, r.ncol))
    idx = deepcopy(r.pointindex)

    # cannot do knn search on empty tree
    if isempty(index)
        return Raster(idx, r.nrow, r.ncol, r.bbox, r.cellsize, r.epsg, r.npoints), trues(r.nrow, r.ncol)
    end

    tree = create_kdtree(cloud; dim=2, index=index)

    if length(tree.data) <= 1
        return r, missings
    end

    @showprogress 1 "Painting missing data..." for icell in 1:length(r)
        isempty(idx[icell]) || continue
        missings[icell] = true
        # for empty cells find nearest neighbor point
        x, y = ind2coords(r.bbox, r.nrow, r.ncol, r.cellsize, icell)
        i, d = knn(tree, [x, y], 1)
        d[1] <= fill_max_dist && push!(idx[icell], index[i[1]])
    end

    # number of unique points unchanged
    Raster(idx, r.nrow, r.ncol, r.bbox, r.cellsize, r.epsg, r.npoints), missings
end


"""filter a raster index with a pointfilter"""
function filter_raster(r::Raster, cloud::Cloud, pointfilter)
    n = length(r)
    # create empty index
    idx = Array{Vector{Int32}}(undef, n) # some stay zero because of skipped points
    for icell in 1:n
        idx[icell] = Int32[]
    end

    # index points that satisfy pointfilter
    npoints = 0
    for icell = 1:n
        for i in r[icell]
            (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
            push!(idx[icell], i)
            npoints += 1
        end
    end

    # return new index
    Raster(idx, r.nrow, r.ncol, r.bbox, r.cellsize, r.epsg, npoints)
end

"generic function to create 2d array from Cloud
inputs
cloud            Cloud
vx, vy        x&y coordinate vectors defining the grid
xy_idx        Cloud indices per gridcell
reducer    statistical function to summarize points per cell
      function reducer(cloud::Cloud, index::Vector{Int})
        Vector{Real}
      end
pointfilter function to filter points
      function pointfilter(cloud::Cloud, index::Int)
        true
      end
"
function rasterize(cloud::Cloud, r::Raster;
    reducer=reducer_minz, # function to calculate value per cell pased on
    pointfilter=nothing, # predicate function to filter individual points
    min_dens=0,
    nodata=-9999.0,
    return_density=false)

    # check dimensions of output
    nlayers = 1 + return_density
    area = r.cellsize^2

    # setup output arrays
    # TODO don't force to Float32
    raster = fill(Float32(nodata), r.nrow, r.ncol, nlayers)

    # loop through cells
    for i in 1:length(r)
        # find points in range
        idx0 = r[i]
        row, col = ind2sub(r, i)
        # loop through points and filter
        if pointfilter != nothing # if filter is given
            idx = Int[] #  indices after filter
            for i in idx0
                pointfilter(cloud, i) || continue
                push!(idx, i)
            end
        else
            idx = idx0 # no filter applied
        end

        np = length(idx)
        # calculate density [points / m2]
        di = Float64(np) / area
        if return_density
            raster[row, col, end] = di # save to end layer
        end
        if min_dens > 0
            di < min_dens && continue # min density threshold
        end

        np < 1 && continue # stat functions don't work on empty arrays
        # set statistics to grid
        raster[row, col, 1:end - return_density] .= reducer(cloud, idx)
    end

    # return raster with statistics, density is saved to last layer
    raster
        end
