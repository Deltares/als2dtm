# runs only using Julia 0.7 at the moment (RV 24 Dec 2020)
# Basic implementation of block percentile tif to tif function.
# Assumes the output resolution divides the input resolution.
# Currently only computes the median, but this is easily changed into percentiles.

using Compat
import Compat.String
using PeatUtils
using ArgParse
using DelimitedFiles
using Statistics

function parse_commandline()
    s = ArgParseSettings(description = "Block percentile raster resampler.")

    @add_arg_table! s begin
        "--res-in"
            help = "input resolution"
            arg_type = Float64
        "--res-out"
            help = "output resolution"
            arg_type = Float64
        "--outdir"
            help = "output directory"
        "--mindatafraction"
            help = "fraction of input data cells required, else nodata"
            arg_type = Float64
            default = 1.0
            range_tester = x -> 0.0 <= x <= 1.0
        "--filelist"
            help = "file with list of lasfiles to process"
            arg_type = String
        "--files"
            help = "list of lasfiles to process"
            arg_type = String
            nargs = 'R'
    end

    return parse_args(s)
end

arg = parse_commandline()

outdir = arg["outdir"]
mindatafraction = arg["mindatafraction"]
res_in = arg["res-in"]
res_out = arg["res-out"]
files = arg["files"]
filelist = arg["filelist"] == nothing ? [] : vec(readdlm(arg["filelist"], '\n', String))
append!(files, filelist)

nsample = round(Int, res_out^2 / res_in^2)
nmin = max(round(Int, nsample * mindatafraction), 1) # at least 1 cell
isdir(outdir) || mkpath(outdir) # ensure it exists

# for block indexing
bl = Int(res_out ÷ res_in) - 1
d = Int(res_out ÷ res_in)

for file in files
    @show file

    rin, x_min, y_max, cellsize, nodata, transform, epsg = read_raster(file)

    nrow1, ncol1 = size(rin)
    # determine new shape: truncated currently
    nrow2, ncol2 = (nrow1 ÷ d, ncol1 ÷ d)
    # create arrays of the same type as the input array
    # initialized to consist fully of nodata values
    rmed = fill(convert(eltype(rin), nodata), nrow2, ncol2)
    rmin = fill(convert(eltype(rin), nodata), nrow2, ncol2)

    isdata(x) = (x != nodata) && !isnan(x)

    for i = 1:nrow2, j = 1:ncol2
        tile = rin[(i-1)*d+1:(i-1)*d+1+bl, (j-1)*d+1:(j-1)*d+1+bl]
        ndata = count(isdata, tile)
        if ndata >= nmin
            pool = filter(isdata, tile)
            rmed[i,j] = median(pool)
            rmin[i,j] = minimum(pool)
        end
    end

    bname = splitext(basename(file))[1]
    tifmed = joinpath(outdir, bname * "-med.tif")
    tifmin = joinpath(outdir, bname * "-min.tif")
    tifdif = joinpath(outdir, bname * "-dif.tif")

    # generate a difference map, preserving nodata
    # both rasters have the same nodata mask
    isnodata = rmin .== nodata
    rdif = rmed - rmin
    rdif[isnodata] .= nodata

    write_raster(tifmed, rmed, x_min, y_max, res_out; epsg=epsg, nodata=-9999)
    write_raster(tifmin, rmin, x_min, y_max, res_out; epsg=epsg, nodata=-9999)
    write_raster(tifdif, rdif, x_min, y_max, res_out; epsg=epsg, nodata=-9999)
end
