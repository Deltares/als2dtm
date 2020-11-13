
using GridOperations
using PeatUtils
using XYZ
using Compat

const overlap = 0.  # width of overlap [m]

# general settings
const epsg = 32748  # UTM 48S
const high_res = 1.0  # cellsize for water grids & and high res dtm [m]
const low_res = 25.0  # cellsize for water grids & and high res dtm [m]
const nodata = -9999.0
const min_points = 0.5

function xyvec(cloud::XYZ.Cloud, sampleratio::Float64)
    n = length(cloud.positions)
    step = round(Int, cld(1, sampleratio))
    # have minimum of points needed
    if n < 10000 | fld(n, step) < 1000
        step = 1
    end
    idx_samples = 1:step:n
    nsample = length(idx_samples)
    A = Vector{Vector{Float64}}() # empy vector of vectors
    sizehint!(A, nsample)
    for i in idx_samples
        p = cloud.positions[i]
        push!(A, [p[1], p[2]])
    end
    A
end

# simplify raster from cloud
function cloud2arr(cloud::XYZ.Cloud, tile_bbox::XYZ.BoundingBox, cellsize::Float64, epsg::Int; pointfilter=nothing, reducer=XYZ.reducer_medz, interpolate=true)
    index = XYZ.define_raster(cloud, tile_bbox, overlap, cellsize; snapgrid=25.0, epsg=epsg, pointfilter=pointfilter)
    raster = XYZ.rasterize(cloud, index; reducer=reducer, min_dens=0)[:,:,1]
    if interpolate
        raster = GridOperations.interp_missing!(copy(raster), eltype(raster)(nodata); mode="kriging")
    end
    raster, index
end

# simplify tif creation
function cloud2tif(folder::String, fn::String, cloud::XYZ.Cloud, tile_bbox::XYZ.BoundingBox, cellsize::Float64, epsg::Int, filter, reducer; interpolate=false)
    raster, index = cloud2arr(cloud, tile_bbox, cellsize, epsg, pointfilter=filter, reducer=reducer, interpolate=interpolate)
    XYZ.grid2tif(folder, fn, index, raster; nodata=nodata)
end


"""pipeline to classify cloud and create dtm"""
function lidar_pipeline(
        filename::String,
        epsg::Int,
        out_dir::String,
        gf_radius::Float64,
        gf_slope::Float64,
        gf_dh_max::Float64,
        gf_dh_min::Float64,
        high_res::Float64,
        low_res::Float64,
        )

    # file admin
    fileid = basename(splitext(filename)[1])
    in_dir = dirname(filename)
    filen = fileid

    out_dir = joinpath(out_dir, fileid)
    isdir(out_dir) || mkpath(out_dir)
    @info("Ouput directory is $(out_dir)")

    # Set up Logging
    @info("Processing $(filen)")
    @info("with ground parameters r:$(gf_radius) \t s:$(gf_slope) \t min:$(gf_dh_min) \t max:$(gf_dh_max)")

    # read las
    @info("--> read $(fileid)")
    @info("--> write results to $(out_dir)")
    cloud, header = XYZ.read_pointcloud(filename)
    @info("  $(cloud)")

    # reset classification
    XYZ.reset_class!(cloud)

    tile_bbox = XYZ.BoundingBox(xmin=header.x_min, ymin=header.y_min,
        xmax=header.x_max, ymax=header.y_max)

    ## CLASSIFY OUTLIERS
    @info("--> 1.) classify outliers")
    # classify low and high outliers
    # evaluate points per in <cellsize> [m] cells
    # low outliers if a jump larger than <dz> [m] in lowest n <max_outliers> [-] points
    # high outliers if points more than <max_height> [m] above lowest non-outlier point
    XYZ.classify_outliers!(cloud;
        cellsize=100.0,   # cellsize used for evaluation of points
        dz=1.0,           # threshold in vertical distance between low points
        max_outliers=15,  # max number of outliers in one cell = low  points to be evaluated
        max_height=100.0)   # threshold in vertical distance for high points

    ## CLASSIFY GROUND with Zhang algorithm. use only non-water, non-outlier, last return points
    @info("--> 1.5) intermediate rasters")

    # make non-filled raster index for high resolution cells with non-outliers and last return points
    r = XYZ.define_raster(cloud, tile_bbox, overlap, high_res; pointfilter=XYZ.unclassified_lastreturn, epsg=epsg)
    if !(r.nrow >= 2 && r.ncol >= 2)
        @info("Defined raster is too small to be classified, skipping")
        @show r
        return
    end

    boundarymask_r = trues(r.nrow, r.ncol)
    smallboundarymask_r = trues(r.nrow, r.ncol)

    zmin = XYZ.rasterize(cloud, r; reducer=XYZ.reducer_minz,  pointfilter=nothing)[:,:,1]
    XYZ.grid2tif(out_dir, "$(filen)_zmin.tif", r, zmin; nodata=nodata)

    # create las mask and dropouts & vegetation binary images
    haspoints = Array{Bool}(XYZ.rasterize(cloud, r; reducer=XYZ.reducer_count)[:,:,1] .> 0)
    las_mask = GridOperations.create_mask(haspoints)
    dropouts = Array{Bool}(las_mask .& broadcast(~, haspoints))  # nodata cells within las mask

    # zmin interpolated surface
    r_filled = XYZ.inpaint_missings_nn(r, cloud)[1]
    zmin_filled = XYZ.rasterize(cloud, r_filled; reducer=XYZ.reducer_minz,  pointfilter=nothing)[:,:,1]
    XYZ.grid2tif(out_dir, "$(filen)_zmin_nn.tif", r, zmin_filled; nodata=nodata)

    ## CLASSIFY GROUND with Zhang algorithm. use only non-water, non-outlier, last return points
    @info("--> 2.) classify terrain (pmf)")

    # apply zhang to grid
    zmax_pmf, flags = GridOperations.pmf_filter(zmin_filled, boundarymask_r, gf_radius, gf_slope, gf_dh_max, gf_dh_min, r.cellsize)
    XYZ.grid2tif(out_dir, "$(filen)_zmax_pmf.tif", r, zmax_pmf; nodata=nodata)
    bin_vegetation = Array{Bool}(flags .> 0) .& broadcast(~, dropouts)

    # use maxz grid for classification
    XYZ.classify_below_surface!(cloud, r, zmax_pmf, 2)
    XYZ.classify_mask!(cloud, r, smallboundarymask_r, 28)  # unused class for edge effects

    # filter zmin grid
    zmin_filtered = copy(zmin)
    zmin_filtered[((bin_vegetation .| dropouts) .| .!smallboundarymask_r)] .= nodata
    XYZ.grid2tif(out_dir, "$(filen)_zmin_filtered.tif", r, zmin_filtered; nodata=nodata)

    # Apply classification to pointcloud
    r_pmf = XYZ.filter_raster(r, cloud, XYZ.ground)
    cloud_min = XYZ.reduce_min(cloud, r_pmf)

    let r_terrain, dtm_lowres, cs = low_res
        r_terrain = XYZ.define_raster(cloud_min, tile_bbox, overlap, cs; snapgrid=low_res, epsg=epsg,
            pointfilter=XYZ.ground)
        dtm_lowres = XYZ.rasterize(cloud_min, r_terrain; reducer=XYZ.reducer_medz,
            min_dens=min_points / (cs^2))[:,:,1]
        XYZ.grid2tif(out_dir, "$(filen)_dtm_$(cs)m.tif", r_terrain, dtm_lowres; nodata=nodata)
    end

    nothing
end
