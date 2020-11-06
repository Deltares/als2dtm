"Decompose a ColorTypes.RGB{N0f16} into red, green and blue color channels"
rgb_components(c::RGB{N0f16}) = red(c), green(c), blue(c)

"write tiffile from Cloud"
function to_tif(
    outDir::String,
    filenames::Vector{String},
    cloud::Cloud,
    cellsize::Real,
    epsg::Int;
    reducer=minz,
    pointfilter=nothing,
    min_dens=0,
    nodata=-9999,
    return_density=false)

    # create grid definition
    r = define_raster(cloud, cellsize; epsg=epsg, pointfilter=pointfilter)

    to_tif(outDir, filenames, cloud, r;
        reducer=reducer, pointfilter=nothing, # filter in define_raster cellsize
        min_dens=min_dens, nodata=nodata, return_density=return_density)
end

"write tiffile from Cloud with given raster definition"
function to_tif(outDir::String, filenames::Vector{String},
    cloud::Cloud, r::Raster;
    reducer=minz, pointfilter=nothing,
    min_dens=0, nodata=-9999, return_density=false)

    # make raster from Cloud based on statistical summary per griddcell and point filter
    grid = rasterize(cloud, r;
        reducer=reducer,
        pointfilter=pointfilter,
        min_dens=min_dens,
        nodata=nodata,
        return_density=return_density)

    # write to file
    grid2tif(outDir, filenames, r, grid;
        nodata=nodata)
end

"""util function to save a grid to tif using properties from raster defintion r
for 2d and 3d (multiple layers) grids; 3d grids are written to multiple files"""
function grid2tif(outDir::String, filenames::Vector{String}, r::Raster, grid::Array; nodata=-9999)
    ndim = ndims(grid)
    grid = copy(grid)
    if ndim == 3
        nlayers = size(grid, 3)
        for i in 1:nlayers
            fn = joinpath(outDir, filenames[i])
            grid[:, :, i][r.mask] = nodata
            write_raster(fn, grid[:, :, i], r.bbox.xmin, r.bbox.ymax, r.cellsize; epsg=XYZ.epsg(r), nodata=nodata)
        end
    elseif ndim == 2
        fn = joinpath(outDir, filenames[1])
        grid[r.mask] = nodata
        write_raster(fn, grid, r.bbox.xmin, r.bbox.ymax, r.cellsize; epsg=XYZ.epsg(r), nodata=nodata)
    end
end
"for 2d arrays"
function grid2tif(outDir::String, filename::String, r::Raster, grid::Array; nodata=-9999)
    @assert (ndims(grid) == 2) "Invalid grid size, use a vector of filenames to save layers to multiple files"
    fn = joinpath(outDir, filename)
    grid = copy(grid)
    grid[r.mask] .= nodata
    write_raster(fn, grid, r.bbox.xmin, r.bbox.ymax, r.cellsize; epsg=XYZ.epsg(r), nodata=nodata)
end

"merge two headers, rejecting incompatible headers, and favoring the first"
function Base.merge(h1::LasHeader, h2::LasHeader)
    # check compatibility
    msg = "Cannot merge files, "
    h1.version_major === h2.version_major || @error(msg * "LAS major version mismatch")
    h1.version_minor === h2.version_minor || @error(msg * "LAS minor version mismatch")
    h1.data_format_id === h2.data_format_id || @error(msg * "Point format mismatch")
    h1.data_format_id === h2.data_format_id || @error(msg * "Point format mismatch")
    h1.data_record_length === h2.data_record_length || @error(msg * "Point record length mismatch")

    # merge selected fields
    records_count = h1.records_count + h2.records_count
    point_return_count = h1.point_return_count + h2.point_return_count
    x_max = max(h1.x_max, h2.x_max)
    x_min = min(h1.x_min, h2.x_min)
    y_max = max(h1.y_max, h2.y_max)
    y_min = min(h1.y_min, h2.y_min)
    z_max = max(h1.z_max, h2.z_max)
    z_min = min(h1.z_min, h2.z_min)

    # Note that because this function is used for appending LAS files, the size of the header
    # must keep the same size. Therefore the VLRs are currently not extended, and neither are
    # the user defined bytes.
    # Furthermore it is assumed that the new coordinates will fit in with the existing h1 scale
    # and offset values.

    # create merged header
    LasHeader(
        h1.file_source_id,
        h1.global_encoding,
        h1.guid_1,
        h1.guid_2,
        h1.guid_3,
        h1.guid_4,
        h1.version_major,
        h1.version_minor,
        h1.system_id,
        h1.software_id,
        h1.creation_doy,
        h1.creation_year,
        h1.header_size,
        h1.data_offset,
        h1.n_vlr,
        h1.data_format_id,
        h1.data_record_length,
        records_count,
        point_return_count,
        h1.x_scale,
        h1.y_scale,
        h1.z_scale,
        h1.x_offset,
        h1.y_offset,
        h1.z_offset,
        x_max,
        x_min,
        y_max,
        y_min,
        z_max,
        z_min,
        h1.variable_length_records, # VLRs currently not extended, would need to be merged but not duplicated
        h1.user_defined_bytes
    )
end

"write LAS points"
function writepoints(io::IO, header::LasHeader, cloud::Cloud)
    n = length(positions(cloud))
    pointtype = LasIO.pointformat(header)
    for i in 1:n
        lasp = laspoint(pointtype, cloud, i, header)
        write(io, lasp)
    end
end

"Construct a LasIO Vector{LasPoint} from a XYZ.Cloud"
function lasio_vector(header::LasHeader, cloud::Cloud)
    n = length(positions(cloud))
    pointtype = LasIO.pointformat(header)
    pointdata = Vector{pointtype}(n)
    for i in 1:n
        pointdata[i] = laspoint(pointtype, cloud, i, header)
    end
    pointdata
end

"write LAS/LAZ file from Cloud"
function to_las(outDir::String, filename::String,
                cloud::Cloud, header::LasHeader)

    fn = joinpath(outDir, filename)

    header = deepcopy(header) # prevent unwanted side effects
    LasIO.update!(header, cloud)

    # converts to FileIO File LAS or LAZ
    # such that it dispatches to the right save
    fio = query(fn)
    save(fio, header, cloud)
    nothing
end
    
# Add save on Cloud be able to write LAS from Cloud without copies
function save(f::File{format"LAS"}, header::LasHeader, cloud::Cloud)
    # convert File to String to get a normal IO object
    # that the other functions can use
    fn = filename(f)
    open(fn, "w") do io
        write(io, "LASF")
        write(io, header)
        writepoints(io, header, cloud)
    end
end

# Add save on Cloud be able to write LAZ from Cloud without copies
function save(f::FileIO.File{FileIO.DataFormat{:LAZ}}, header::LasHeader, cloud::Cloud)
    # pipes las to laszip to write laz
    fn = filename(f)
    open(`laszip -olaz -stdin -o "$fn"`, "w") do s
        savebuf(s, header, cloud)
    end
end

# Add savebuf on Cloud be able to write LAZ from Cloud without copies
function savebuf(s::IO, header::LasHeader, cloud::Cloud)
    # checks
    header_n = header.records_count
    n = length(cloud.positions)
    msg = "number of records in header ($header_n) does not match data length ($n)"
    @assert header_n == n msg
    pointtype = LasIO.pointformat(header)

    # write header
    write(s, magic(format"LAS"))
    write(s, header)

    # 2048 points seemed to be an optimum for the libLAS_1.2.las testfile
    npoints_buffered = 2048
    bufsize = header.data_record_length * npoints_buffered
    buf = IOBuffer(bufsize)
    # write points
    for i in 1:n
        p = laspoint(pointtype, cloud, i, header)
        write(buf, p)
        if rem(i, npoints_buffered) == 0
            write(s, take!(buf))
        elseif i == n
            write(s, take!(buf))
        end
    end
            end

"general XYZ (csv) writer"
function to_xyz(outDir::String, filename::String, cloud::Cloud;
    precision=2, delimiter=",",
    attributes=Symbol[],
    pointfilter=nothing)

    fn = joinpath(outDir, filename)
    nrow = length(positions(cloud))
    write_attributes = length(attributes) > 0
    header =  xyz_header(attributes; delimiter=delimiter)

    # open file
    open(fn, "w") do f
        # write header
        write(f, header)

        # loop through points and write a line per point to file
        for i in 1:nrow
            (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
            line = String[]
            # read and position point data and format
            for p in positions(cloud)[i]
                push!(line, format(p, precision=precision))
                push!(line, delimiter)
            end
            # read and attribute point data and format
            if write_attributes
                for key in attributes
                    if typeof(cloud[key][i]) <: Integer
                        push!(line, format(Int(cloud[key][i]), precision=0))
                    else
                        push!(line, format(cloud[key][i], precision=precision))
                    end
                    push!(line, delimiter)
                end
            end

          # write formatted line
          write(f, string(join(line)[1:end - length(delimiter)], "\n"))
        end
    end # close file
end

"XYZ (csv) writer for points along profile. The coordinates in profile orientation
are written in extra columns px and py to xyz file"
function to_xyz(outDir::String, filename::String, cloud::Cloud, prof::Profile;
    precision=2, delimiter=",",
    attributes=Symbol[],
    pointfilter=nothing)

    fn = joinpath(outDir, filename)
    nrow = length(prof)
    write_attributes = length(attributes) > 0
    header = []
    header =  xyz_header(cat(1, ["px", "py"], attributes); delimiter=delimiter)

    # open file
    open(fn, "w") do f
        # write header
        write(f, header)

        # loop through points and write a line per point to file
        ip = 0
        for i in pointindex(prof)
            ip += 1
            (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
            line = String[] # create vector with formated point data
            # read profile coordinates and format
            # read and position point data and format
            for p in positions(cloud)[i]
                push!(line, format(p, precision=precision))
                push!(line, delimiter)
            end
                pos = Float64[pointlx(prof)[ip], pointly(prof)[ip]]
            for p in pos
                push!(line, format(p, precision=precision))
                push!(line, delimiter)
            end
                # read and attribute point data and format
            if write_attributes
                for key in attributes
                    if typeof(cloud[key][i]) <: Integer
                        push!(line, format(Int(cloud[key][i]), precision=0))
                    else
                        push!(line, format(cloud[key][i], precision=precision))
                    end
                  push!(line, delimiter)
                end
            end

            # write formatted line
            write(f, string(join(line)[1:end - length(delimiter)], "\n"))
        end
    end # close file
end

"create header for xyz files"
function xyz_header(attributes::Vector{}; delimiter=",")
    header = ["x,","y,","z,"]
    if length(attributes) > 0
        for key in attributes
            push!(header, string(key))
            push!(header, delimiter)
        end
    end
    string(join(header)[1:end - length(delimiter)], "\n")
end
