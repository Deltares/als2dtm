using XYZ
using StaticArrays
using Base.Test

# use test file from LasIO
filename = "libLAS_1.2.las" # point format 0
filepath = joinpath(Pkg.dir("LasIO"), "test", filename)
workdir = dirname(@__FILE__)

@testset "Cloud contruction and indexing" begin
    cloud, header = XYZ.read_pointcloud(filepath)

    attr = cloud.attributes
    @test isa(attr, Dict{Symbol,Vector})
    @test Set(keys(attr)) == Set(Symbol[:number_of_returns, :edge_of_flight_line, :intensity,
        :key_point, :scan_direction, :withheld, :pt_src_id, :classification, :scan_angle, :synthetic,
        :return_number, :user_data])

    @test XYZ.boundingbox(cloud) == XYZ.calc_bbox(cloud)

    cloudp1 = cloud[1:1]
    @test length(cloudp1) == 1
    @test XYZ.positions(cloudp1) == [SVector{3,Float64}(1.44013394e6, 375000.23, 846.66)]
    @test cloudp1[:intensity] == [0x00fa]
    @test cloudp1[:scan_angle] == [Int8(0)]
    @test cloudp1[:user_data] == [0x00]
    @test cloudp1[:pt_src_id] == [0x001d]
    @test cloudp1[:return_number] == [0x00]
    @test cloudp1[:number_of_returns] == [0x00]
    @test cloudp1[:scan_direction] == [false]
    @test cloudp1[:edge_of_flight_line] == [false]
    @test cloudp1[:classification] == [0x02]
    @test cloudp1[:synthetic] == [false]
    @test cloudp1[:key_point] == [false]
    @test cloudp1[:withheld] == [false]
end

@testset "Outlier classification" begin
    cloud, header = XYZ.read_pointcloud(filepath)

    XYZ.reset_class!(cloud)
    @test all(c -> c == 0x00, cloud[:classification])

    XYZ.classify_outliers!(cloud;
        cellsize=100.0,   # cellsize used for evaluation of points
        dz=1.0,           # threshold in vertical distance between low points
        max_outliers=15,  # max number of outliers in one cell = low  points to be evaluated
        max_height=100.0) # threshold in vertical distance for high points

    low_noise = UInt8(7) # from ASPRS Standard LIDAR Point Classes
    high_noise = UInt8(18) # from ASPRS Standard LIDAR Point Classes (LAS 1.4)
    @test count(c -> c == low_noise, cloud[:classification]) == 13
    @test count(c -> c == high_noise, cloud[:classification]) == 1948
end

@testset "Rasterization" begin
    cloud, header = XYZ.read_pointcloud(filepath)
    cellsize = 10.0
    nground = 402812
    nodata = -9999.0
    rasterindex = XYZ.define_raster(cloud, cellsize)
    @test rasterindex.npoints == 497536
    @test typeof(rasterindex) == XYZ.Raster{Int32}
    @test rasterindex.nrow == rasterindex.ncol == 500

    raster = XYZ.rasterize(cloud, rasterindex; reducer=XYZ.reducer_minz)
    @test size(raster) == (500, 500, 1)
    @test eltype(raster) == Float32
    @test count(x -> x == nodata, raster) == 3336

    rasterindex_ground = XYZ.filter_raster(rasterindex, cloud, XYZ.ground)
    @test rasterindex_ground.npoints == nground

    # pointfilter
    @test count(i -> XYZ.ground(cloud, i), 1:length(cloud)) == nground
end
