__precompile__()

module XYZ

# export
#     # cloud
#     Cloud, BoundingBox, positions, attributes, boundingbox, calc_bbox, update_bbox!,
#     # pointfilters
#     ground, water, notground, notoutlier, lastreturn,
#     # classifiers
#     classify_ground!, classify_ground_approx!, classify_outliers!, classify_water!,
#     classify_min!, classify_max!, classify!, reset_class!,
#     classify_below_surface!, classify_above_surface!, classify_water_approx!,
#     # raster
#     Raster, rasterize, define_raster,
#     # profile
#     define_profile, filter_profile,
#     # reducers
#     reduce_min, reduce_max, reduce_pointfilter,
#     # writers
#     to_tif, to_las, to_xyz, grid2tif,
#     # readers
#     read_pointcloud,
#     # utils
#     describe, paint!, create_kdtree, point_density!, grid2cloud_attribute!, copy_classification!

using Compat
using Nullables
using ColorTypes
using FixedPointNumbers
using Colors
using GDAL
using FileIO
using LasIO
using NearestNeighbors
using PeatUtils
using ProgressMeter
using SortingAlgorithms
using StaticArrays

include("cloud.jl")
include("pointfilters.jl")
include("rasterize.jl")
include("cloudutils.jl")
include("profiler.jl")
include("cloudfilter.jl")
include("cloudclassifier.jl")
include("lasio.jl")
include("writers.jl")

end
