using als2dtm

using ArgParse

function parse_commandline()
    s = ArgParseSettings(description="Create a DTM from a pointcloud.")

    @add_arg_table! s begin
        "--epsg"
            help = "CRS EPSG code"
            arg_type = Int
            default = 32748
        "--slope"
            help = "PMF slope"
            arg_type = Float64
            default = 0.6
            range_tester = x -> x > 0.0
        "--dh_min"
            help = "PMF dh_min"
            arg_type = Float64
            default = 0.5
            range_tester = x -> x > 0.0
        "--dh_max"
            help = "PMF dh_max"
            arg_type = Float64
            default = 1.1
            range_tester = x -> x > 0.0
        "--high_res"
            help = "High resolution to rasterize"
            arg_type = Float64
            default = 1.0
            range_tester = x -> x > 0.0
        "--low_res"
            help = "Low resolution to rasterize"
            arg_type = Float64
            default = 100.
            range_tester = x -> x > 0.0
        "--radius"
            help = "PMF radius in coordinate system units"
            arg_type = Float64
            default = 16.
            range_tester = x -> x > 0.0
        "--infile", "-i"
            help = "path to input .las file"
            arg_type = String
            range_tester = x -> isfile(x)
            required = true
        "--outfolder", "-o"
            help = "path to output folder"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

arg = parse_commandline()
epsg = arg["epsg"]
slope = arg["slope"]
dh_min = arg["dh_min"]
dh_max = arg["dh_max"]
high_res = arg["high_res"]
low_res = arg["low_res"]
infile = arg["infile"]
outfolder = arg["outfolder"]
radius = arg["radius"]

isdir(outfolder) || mkpath(outfolder) # ensure it exists

# infile = "data/small.laz"
# infile = "area1.las"
# outfolder = "output/"

# # PMF settings
# radius = 16.
# slope = 0.6
# dhmax = 1.1
# dhmin = 0.5

# # Resolution settings
# high_res = 1.0
# low_res = 100.0

lidar_pipeline(infile, epsg, outfolder, radius, slope, dh_max, dh_min, high_res, low_res)
