# Interpolate missing data in rasters. Currently only Inverse Distance Weighting
# is implemented, using a maxium search radius and it is blocking. This means
# that data cells hidden behind other data cells are ignored.

using IDWInterpolations
using ArgParse
using PeatUtils

function parse_commandline()
    s = ArgParseSettings(description="Interpolate missing data in rasters.")

    @add_arg_table! s begin
        "--nodata"
            help = "nodata value"
            arg_type = Float64
            default = -9999.0
        "--power"
            help = "power used for distance weighting"
            arg_type = Float64
            default = 2.0
            range_tester = x -> x >= 1.0
        "--radius"
            help = "maximum search distance for data points, in projection units"
            arg_type = Float64
            default = 10_000.0 # 10 km assuming meters
            range_tester = x -> x > 0.0
        "--infile", "-i"
            help = "path to input raster"
            arg_type = String
            range_tester = x -> isfile(x)
            required = true
        "--outfile", "-o"
            help = "path to output raster"
            arg_type = String
            required = true
    end

return parse_args(s)
end

arg = parse_commandline()
nodata = arg["nodata"]
power = arg["power"]
fname_in = arg["infile"]
fname_out = arg["outfile"]
searchradius_m = arg["radius"]

outdir = dirname(fname_out)
isdir(outdir) || mkpath(outdir) # ensure it exists
isfile(fname_out) && rm(fname_out) # delete existing output files

rd, x_min, y_max, cellsize, nodata_in, transform, epsg = read_raster(fname_in)

# change any NaN to the given nodata value (easier with PCRaster)
rd[isnan.(rd)] .= nodata
# println(rd)
println(size(rd))
println(nodata)
println(nodata_in)
# the supplied nodata value overrides the input raster nodata value
rd[rd .== nodata_in] .= nodata
searchradius_cells = searchradius_m / cellsize

ai = blocking_idw(rd, nodata, searchradius_cells, power)

write_raster(fname_out, ai, x_min, y_max, cellsize; epsg=epsg, nodata=nodata)
