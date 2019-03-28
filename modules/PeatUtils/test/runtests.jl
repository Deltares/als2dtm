using PeatUtils
using Base.Test
using GDAL

# reuse GDAL.jl test data
tifpath = joinpath(Pkg.dir("GDAL"), "test", "data", "utmsmall.tif")
tifpath_local = joinpath(dirname(@__FILE__), "data", "utmsmall.tif")

cp(tifpath, tifpath_local, remove_destination=true)

@test getnodata(tifpath) == -1.0e10

nodata = -9999.0
setnodata!(tifpath_local, nodata)
@test getnodata(tifpath_local) == nodata

# test pointsample()
let
    GDAL.allregister()
    dataset = GDAL.open(tifpath_local, GDAL.GA_ReadOnly)
    geotransform = zeros(6)
    GDAL.getgeotransform(dataset, geotransform)
    band = GDAL.getrasterband(dataset, 1)

    @test pointsample(band, geotransform, 443334.9, 3750015.9) === UInt8(189)
    @test pointsample(band, geotransform, 443167.9, 3749928.3) === UInt8(173)

    GDAL.close(dataset)
    GDAL.destroydrivermanager()
end


A, x_min, y_max, cellsize, epsg, nodata = read_raster(tifpath_local)

@test x_min === 440720.0
@test y_max === 3.75132e6
@test cellsize === 60.0
@test epsg === Nullable{Int64}(26711)
@test nodata === -9999.0

tifpath_local_out = joinpath(dirname(@__FILE__), "data", "utmsmall_out.tif")
write_raster(tifpath_local_out, A, x_min, y_max, cellsize; epsg=epsg)

ref_wkt_utm48n = """PROJCS["WGS 84 / UTM zone 48N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32648"]]"""
ref_epsg_utm48n = Nullable(32648)
@test epsg2wkt(ref_epsg_utm48n) == ref_wkt_utm48n
@test wkt2epsg(ref_wkt_utm48n) === ref_epsg_utm48n

rm(tifpath_local)
rm(tifpath_local_out)

# test without projection information
# A, x_min, y_max, cellsize, epsg, nodata = read_raster("a.tif")
# @test isnull(epsg)
# @test isa(epsg, Nullable{Int})
