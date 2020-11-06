#= 
TODO
update header based on Clouds =#


"""Prepare and allocate an empty dict to be used as a Cloud attribute table
based on a list of needed attributes"""
function prepare_attributes(pointtype::Type{T}, n::Integer) where T <: LasPoint
    n = Int(n)  # cannot construct BitVector(n) with n <: UInt32 (fix in julia)
    # create dict with all Las point format 0 properties
    attr = Dict{Symbol,Vector}(
        :intensity => zeros(UInt16, n),
        :return_number => zeros(UInt8, n),
        :number_of_returns => zeros(UInt8, n),
        :scan_direction => zeros(Bool, n),
        :edge_of_flight_line => zeros(Bool, n),
        :classification => zeros(UInt8, n),
        :synthetic => zeros(Bool, n),
        :key_point => zeros(Bool, n),
        :withheld => zeros(Bool, n),
        :scan_angle => zeros(Int8, n), # -90 to +90
        :user_data => zeros(UInt8, n),
        :pt_src_id => zeros(UInt16, n),
    )
    # add additional attributes based on availability
    if pointtype <: LasIO.LasPointTime
        attr[:time] = zeros(Float64, n)
    end
    if pointtype <: LasIO.LasPointColor
        attr[:color] = Vector{RGB{N0f16}}(n)
    end
    attr
end

# this function is a performance bottleneck
function add_point!(
    lasp::LasPoint,
    header::LasHeader,
    pointdata::Vector{SVector{3,Float64}},
    attr::Dict{Symbol,Vector},
    i::Integer,
    pointtype::Type{T}) where T <: LasPoint

        pointdata[i] = SVector{3,Float64}(
        xcoord(lasp, header),
        ycoord(lasp, header),
        zcoord(lasp, header))
    attr[:intensity][i] = intensity(lasp)
    attr[:return_number][i] = return_number(lasp)
    attr[:number_of_returns][i] = number_of_returns(lasp)
    attr[:scan_direction][i] = scan_direction(lasp)
    attr[:edge_of_flight_line][i] = edge_of_flight_line(lasp)
    attr[:classification][i] = classification(lasp)
    attr[:synthetic][i] = synthetic(lasp)
    attr[:key_point][i] = key_point(lasp)
    attr[:withheld][i] = withheld(lasp)
    attr[:scan_angle][i] = reinterpret(Int8, scan_angle(lasp))  # fix laszip reader
    attr[:user_data][i] = user_data(lasp)
    attr[:pt_src_id][i] = pt_src_id(lasp)
    # add additional attributes based on availability
    if pointtype <: LasIO.LasPointTime
        attr[:time][i] = gps_time(lasp)
    end
    if pointtype <: LasIO.LasPointColor
        attr[:color][i] = RGB(lasp)
    end
end

function read_pointcloud(s::Union{Stream{format"LAS"},Base.Process})
    LasIO.skiplasf(s)
    header = read(s, LasHeader)

    n = Int(header.records_count)
    pointtype = LasIO.pointformat(header)

    # pre allocate the parts that will go into the Cloud
    pointdata = Vector{SVector{3,Float64}}(undef, n)
    attr = prepare_attributes(pointtype, n)

    @showprogress 1 "Reading pointcloud..." for i = 1:n
        lasp = read(s, pointtype)
        add_point!(lasp, header, pointdata, attr, i, pointtype)
    end
    header, pointdata, attr
end

function read_pointcloud(filepath::File{format"LAS"})
    header, pointdata, attr = open(filepath) do io
        read_pointcloud(io)
    end
    cloud = Cloud(pointdata, attr)
    cloud, header
end

function read_pointcloud(filepath::File{format"LAZ"})
    fn = filename(filepath)
    header, pointdata, attr = open(`laszip -olas -stdout -i "$fn"`) do s
        read_pointcloud(s)
    end
    cloud = Cloud(pointdata, attr)
    cloud, header
end

"""
    read_pointcloud(filepath::AbstractString)

Read a pointcloud stored in LAS or LAZ format to a Cloud
"""
function read_pointcloud(filepath::AbstractString)
    # converts to FileIO File LAS or LAZ
    # such that it dispatches to the right read_pointcloud
    fio = query(filepath)
    read_pointcloud(fio)
end

function laspoint_shared(cloud::Cloud, i::Integer, h::LasHeader)
    # fields shared among all point formats
    pos = positions(cloud)[i]
    x = xcoord(pos[1], h)
    y = ycoord(pos[2], h)
    z = zcoord(pos[3], h)
    intensity = cloud[:intensity][i]
    flagb = flag_byte(cloud[:return_number][i], cloud[:number_of_returns][i],
        cloud[:scan_direction][i], cloud[:edge_of_flight_line][i])
    rawcl = raw_classification(cloud[:classification][i], cloud[:synthetic][i],
        cloud[:key_point][i], cloud[:withheld][i])
    scan_angle = cloud[:scan_angle][i]
    user_data = cloud[:user_data][i]
    pt_src_id = cloud[:pt_src_id][i]

    x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id
end

"Cloud point to LasPoint0"
function laspoint(pointtype::Type{LasPoint0}, cloud::Cloud, i::Integer, h::LasHeader)
    x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id = laspoint_shared(cloud, i, h)
    LasPoint0(x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id)
end

"Cloud point to LasPoint1"
function laspoint(pointtype::Type{LasPoint1}, cloud::Cloud, i::Integer, h::LasHeader)
    x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id = laspoint_shared(cloud, i, h)
    gpst = cloud[:time][i]
    LasPoint1(x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id,
        gpst)
end

"Cloud point to LasPoint2"
function laspoint(pointtype::Type{LasPoint2}, cloud::Cloud, i::Integer, h::LasHeader)
    x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id = laspoint_shared(cloud, i, h)
    r, g, b = rgb_components(cloud[:color][i])
    LasPoint2(x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id,
        r, g, b)
end

"Cloud point to LasPoint3"
function laspoint(pointtype::Type{LasPoint3}, cloud::Cloud, i::Integer, h::LasHeader)
    x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id = laspoint_shared(cloud, i, h)
    gpst = cloud[:time][i]
    r, g, b = rgb_components(cloud[:color][i])
    LasPoint3(x, y, z, intensity, flagb, rawcl, scan_angle, user_data, pt_src_id,
        gpst, r, g, b)
end

"Update the header bounding box and counts based on point data.
Extends the LasIO.update! to also work on `Cloud`."
function LasIO.update!(h::LasHeader, cloud::Cloud)
    n = length(cloud)
    x_min, y_min, z_min = Inf, Inf, Inf
    x_max, y_max, z_max = -Inf, -Inf, -Inf
    point_return_count = zeros(UInt32, 5)
    for (i, p) in enumerate(positions(cloud))
        x, y, z = p
        if x < x_min
            x_min = x
        end
        if y < y_min
            y_min = y
        end
        if z < z_min
            z_min = z
        end
        if x > x_max
            x_max = x
        end
        if y > y_max
            y_max = y
        end
        if z > z_max
            z_max = z
        end
        # add max statemate to ensure no accidental zeros are passed
        point_return_count[max(cloud[:return_number][i], 1)] += 1
    end
    h.x_min = x_min
    h.y_min = y_min
    h.z_min = z_min
    h.x_max = x_max
    h.y_max = y_max
    h.z_max = z_max
    h.records_count = n
    h.point_return_count = point_return_count
    nothing
end
