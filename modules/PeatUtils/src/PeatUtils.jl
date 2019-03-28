__precompile__()

module PeatUtils
using Nullables
using Compat
import Compat.String
using JSON
using SortingAlgorithms
# using DataArrays
using GeoJSON
using GDAL

export getnodata, setnodata!, rio_info
export wkt, gdalarray
export uint_mapping
export fraction_na
export pointsample
export driver_list
export outside
export read_raster
export write_raster
export wkt2epsg, epsg2wkt

"""
Load a `NULL`-terminated list of strings

That is it expects a "StringList", in the sense of the CPL functions, as a
NULL terminated array of strings.

Function copied over from ArchGDAL.jl
"""
function unsafe_loadstringlist(pstringlist::Ptr{Cstring})
    stringlist = Vector{String}()
    (pstringlist == C_NULL) && return stringlist
    i = 1
    item = unsafe_load(pstringlist, i)
    while Ptr{UInt8}(item) != C_NULL
        push!(stringlist, unsafe_string(item))
        i += 1
        item = unsafe_load(pstringlist, i)
    end
    stringlist
end

"Get the WKT of an Integer EPSG code"
function epsg2wkt(epsg::Nullable{Int})
    if isnull(epsg)
        return "" # missing projections are represented as empty strings
    else
        epsgcode = get(epsg)
        srs = GDAL.newspatialreference(C_NULL)
        GDAL.importfromepsg(srs, epsgcode)
        wkt_ptr = Ref(Cstring(C_NULL))
        GDAL.exporttowkt(srs, wkt_ptr)
        return unsafe_string(wkt_ptr[])
    end
end

epsg2wkt(epsg::Integer) = epsg2wkt(Nullable{Int}(epsg))

"""For EPSG strings like "4326" or "EPSG:4326" """
function epsg2wkt(epsg::String)
    if isempty(epsg)
        return ""
    end
    i = findlast(epsg, ':') + 1 # also works if : is not there
    epsgcode = Nullable{Int}(parse(Int, epsg[i:end]))
    epsg2wkt(epsgcode)
end

"Get the Nullable{Int} EPSG code from a WKT string"
function wkt2epsg(wkt::String)
    if isempty(wkt)
        return Nullable{Int}() # no projection
    else
        # no projection
        srs = GDAL.newspatialreference(C_NULL)
        GDAL.importfromwkt(srs, [wkt])
        epsg = parse(Int, GDAL.getauthoritycode(srs, C_NULL))
        return Nullable{Int}(epsg)
    end
end

# "Convert a DataArray to an ordinary array inplace"
# function convert!(Array, da::DataArray, nodata)
#     da.data[da.na] = nodata
#     da.data
# end

"fallthrough"
convert!(Array, arr::AbstractArray, nodata) = arr

# turn 2D DataArray into regular 1 band array
function gdalarray(da::AbstractArray{T, 2}, nodata) where T
    if eltype(da) in (Int32, Int64) # not supported by GDAL
        GDT, totype = GDAL.GDT_Int32, Int32
    elseif eltype(da) == Float32
        GDT, totype = GDAL.GDT_Float32, Float32
    else # fallthrough: convert to Float64
        GDT, totype = GDAL.GDT_Float64, Float64
    end
    arr = convert!(Array, da, nodata)  # see fallthrough above
    assert(eltype(da) == totype)
    arr, GDT
end

const gdt_lookup = Dict{DataType, GDAL.GDALDataType}(
    UInt8 => GDAL.GDT_Byte,
    UInt16 => GDAL.GDT_UInt16,
    Int16 => GDAL.GDT_Int16,
    UInt32 => GDAL.GDT_UInt32,
    Int32 => GDAL.GDT_Int32,
    Float32 => GDAL.GDT_Float32,
    Float64 => GDAL.GDT_Float64
)

const type_lookup = Dict{GDAL.GDALDataType, DataType}(
    GDAL.GDT_Byte => UInt8,
    GDAL.GDT_UInt16 => UInt16,
    GDAL.GDT_Int16 => Int16,
    GDAL.GDT_UInt32 => UInt32,
    GDAL.GDT_Int32 => Int32,
    GDAL.GDT_Float32 => Float32,
    GDAL.GDT_Float64 => Float64
)


"Get the nodata value from the raster metadata"
function getnodata(fname::AbstractString)
    GDAL.allregister()
    dataset = GDAL.open(fname, GDAL.GA_ReadOnly)
    band = GDAL.getrasterband(dataset, 1)
    nodata = GDAL.getrasternodatavalue(band, C_NULL)
    GDAL.close(dataset)
    GDAL.destroydrivermanager()
    nodata
end

"Set a new nodata value in the raster metadata"
function setnodata!(fname::AbstractString, nodata::Real)
  GDAL.allregister()
  dataset = GDAL.open(fname, GDAL.GA_Update)
  band = GDAL.getrasterband(dataset, 1)
  GDAL.setrasternodatavalue(band, nodata)
  GDAL.close(dataset)
  GDAL.destroydrivermanager()
  nodata
end

"Returns the output of rio info as a Dict"
function rio_info(fname::AbstractString)
    infojson = readall(`rio info $fname`)
    JSON.parse(infojson)
end

"SortingAlgorithms.jl only has uint_mapping(o::Perm, i::Int)"
SortingAlgorithms.uint_mapping(o::Base.Order.Perm{Base.Order.ForwardOrdering,Array{UInt32,1}}, i::UInt32) = SortingAlgorithms.uint_mapping(o.order, o.data[i])

# "Returns the fraction of NA cells in a DataArray"
# fraction_na(da::DataArray) = sum(da.na) / length(da.na)

function driver_list()
    driverlist = String[]
    for i = 0:(GDAL.getdrivercount() - 1)
        driver = GDAL.getdriver(i)
        if driver != C_NULL
            push!(driverlist, GDAL.getdrivershortname(driver))
        end
    end
    return driverlist
end

"check if a point falls outside a bounding box"
function outside(x_min, y_min, x_max, y_max, x, y)
    x < x_min || y < y_min || x > x_max || y > y_max
end

"check if a point falls outside a bounding box"
function outside(bbox, x, y)
    x < bbox[1] || y < bbox[2] || x > bbox[3] || y > bbox[4]
end

"Indexing with points into a GDAL RasterBand, use like this: band[x, y]"
function pointsample(band::Ptr{GDAL.GDALRasterBandH}, geotransform::Vector{Float64}, x::AbstractFloat, y::AbstractFloat)
    x_min = geotransform[1]
    x_res = geotransform[2]
    y_max = geotransform[4]
    y_res = geotransform[6]

    gdt = GDAL.getrasterdatatype(band)
    T = type_lookup[gdt]
    valref = Ref(zero(T)) # black

    nxoff = Int(fld(x - x_min, x_res))
    nyoff = Int(fld(y - y_max, y_res))

    GDAL.rasterio(band, GDAL.GF_Read, nxoff, nyoff, 1, 1,
                  valref, 1, 1, gdt, 0, 0)
    valref[]
end

"Indexing with points into a GDAL RasterBand, use like this: band[x, y]"
function pointsample(T::DataType, band::Ptr{GDAL.GDALRasterBandH}, geotransform::Vector{Float64}, x::AbstractFloat, y::AbstractFloat)
    x_min = geotransform[1]
    x_res = geotransform[2]
    y_max = geotransform[4]
    y_res = geotransform[6]

    gdt = gdt_lookup[T]
    valref = Ref(zero(T)) # black

    nxoff = Int(fld(x - x_min, x_res))
    nyoff = Int(fld(y - y_max, y_res))

    GDAL.rasterio(band, GDAL.GF_Read, nxoff, nyoff, 1, 1,
                  valref, 1, 1, gdt, 0, 0)
    valref[]
end

"Get the data, origin"
function read_raster(fname::String)
    GDAL.allregister()
    dataset = GDAL.open(fname, GDAL.GA_ReadOnly)
    band = GDAL.getrasterband(dataset, 1) # supports single band only
    nodata = GDAL.getrasternodatavalue(band, C_NULL)

    # initialize array to read in with correct dimensions and datatype
    xsize = GDAL.getrasterxsize(dataset)
    ysize = GDAL.getrasterysize(dataset)
    gdt = GDAL.getrasterdatatype(band)
    arrtype = type_lookup[gdt]
    A = zeros(arrtype, xsize, ysize) # dimensions reversed such that we can transpose it back
    # read complete band
    GDAL.rasterio(band, GDAL.GF_Read, 0, 0, size(A,1), size(A,2),
        A, size(A,1), size(A,2), gdt, 0, 0)

    transform = zeros(6)
    GDAL.getgeotransform(dataset, transform)
    x_min = transform[1]
    y_max = transform[4]
    cellsize = transform[2]

    GDAL.close(dataset)
    GDAL.destroydrivermanager()

    Float64.(A'), x_min, y_max, cellsize, nodata, transform
end

"Create a GDAL raster dataset"
function create_raster(
    fname::AbstractString,
    A::Matrix{T},
    x_min::Real,
    y_max::Real,
    cellsize::Real,
    gdaldriver::Ptr{GDAL.GDALDriverH};
    epsg::Nullable{Int}=Nullable{Int}(),
    nodata::Union{Real, Void}=nothing) where T

    bandcount = 1 # this function supports only 1 band rasters
    gdt = gdt_lookup[T]

    # Set compression options for GeoTIFFs
    if gdaldriver == GDAL.getdriverbyname("GTiff")
        options = ["COMPRESS=DEFLATE","TILED=YES"]
    else
        options = String[]
    end

    dstdataset = GDAL.create(gdaldriver, fname, size(A,2), size(A,1), bandcount, gdt, options)

    transform = Float64[x_min, cellsize, 0.0, y_max, 0.0, -cellsize]
    GDAL.setgeotransform(dstdataset,transform)

    projection = epsg2wkt(epsg)
    GDAL.setprojection(dstdataset, projection)

    dstband = GDAL.getrasterband(dstdataset,1)
    GDAL.rasterio(dstband, GDAL.GF_Write, 0, 0, size(A,2), size(A,1),
        A', size(A,2), size(A,1), gdt, 0, 0)

    nodata == nothing || GDAL.setrasternodatavalue(dstband, nodata)
    dstdataset
end

"Write an array to a geospatial GDAL raster"
function write_raster(fname::String, A::Matrix{T}, x_min::Real, y_max::Real, cellsize::Real;
                         epsg::Nullable{Int}=Nullable{Int}(),
                         driver::String="GTiff",
                         nodata::Union{Real, Void}=nothing) where T

     # trying to catch some errors before entering unsafe GDAL territory
     # for safer alternatives look at ArchGDAL.jl
     dir = dirname(fname)
     if isfile(fname)
         rm(fname)
     elseif !isdir(dir)
         mkpath(dir)
     end

    GDAL.allregister()
    gdaldriver = GDAL.getdriverbyname(driver)

    if GDAL.getmetadataitem(gdaldriver, "DCAP_CREATE", "") == "YES"
        ds = create_raster(fname, A, x_min, y_max, cellsize, gdaldriver; epsg=epsg, nodata=nodata)
        GDAL.close(ds)
    elseif GDAL.getmetadataitem(driver, "DCAP_CREATECOPY", "") == "YES"
        memdriver = GDAL.getdriverbyname("MEM")
        dsmem = create_raster(fname, A, x_min, y_max, cellsize, memdriver; epsg=epsg, nodata=nodata)
        progressfunc = convert(Ptr{GDAL.GDALProgressFunc}, C_NULL)
        ds = GDAL.createcopy(gdaldriver, fname, dsmem, 0, C_NULL, progressfunc, C_NULL)
        GDAL.close(dsmem)
        GDAL.close(ds)
    else
        throw(DomainError("GDAL driver $driver does not support CREATE or CREATECOPY"))
    end

    GDAL.destroydrivermanager()
    nothing
end

"Write an array to a geospatial GDAL raster"
function write_raster(fname::String, A::Matrix{T}, x_min::Real, y_max::Real,
                         cellsize::Real, epsg::Int;
                         kwargs...) where T
    write_raster(fname, A, x_min, y_max, cellsize; epsg=Nullable{Int}(epsg), kwargs...)
end

# function write_raster(fname::String, A::DataMatrix{T}, x_min::Real, y_max::Real, cellsize::Real;
#                          nodata::Union{Real, Void}=nothing,
#                          kwargs...) where T
#      ndval = nodata == nothing ? 0 : nodata
#      # convert missing data to a normal Matrix
#      B = convert(Matrix, A, ndval)
#      write_raster(fname, B, x_min, y_max, cellsize; nodata=ndval, kwargs...)
# end

end
