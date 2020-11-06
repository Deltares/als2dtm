__precompile__()

module IDWInterpolations

using NearestNeighbors
using Distances
using ProgressMeter
using LinearAlgebra

export blocking_idw


function to_xy_z(arr, nodata)
    n = length(arr)
    nrow, ncol = size(arr)
    # determine number of data cells
    ndata = sum(arr .!= nodata)
    xy = zeros(2, ndata)
    z = zeros(ndata)
    ind = 1
    for i in 1:nrow, j in 1:ncol
        arr[i, j] == nodata && continue
        xy[:, ind] = Float64[i, j]
        z[ind] = arr[i, j]
        ind += 1
    end
    xy, z
end

function isbarricaded(barriers, x1::Int, y1::Int, x2::Int, y2::Int)
    dx = x2 - x1
    dy = y2 - y1

    # rotate line if it is steep
    is_steep = abs(dy) > abs(dx)
    if is_steep
        x1, y1 = y1, x1
        x2, y2 = y2, x2
    end

    # swap start and end points if necessary and store swap state
    swapped = false
    if x1 > x2
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = true
    end

    # recalculate differentials
    dx = x2 - x1
    dy = y2 - y1

    # calculate error
    error = trunc(Int, dx / 2.0)
    ystep = ifelse(y1 < y2, 1, -1)

    # iterate over bounding box generating points between start and end
    y = y1
    for x in x1:x2
        if is_steep
            coord1, coord2 = y, x
        else
            coord1, coord2 = x, y
        end
        # if it is a barrier, but not the beginning or end point
        if barriers[coord1, coord2] && !([x, y] == [x1, y1]) && !([x, y] == [x2, y2])
            return true
        end
        error -= abs(dy)
        if signbit(error)
            y += ystep
            error += dx
        end
    end
    false
end

# calculate the Inverse Distance Weighting interpolation for one pixel
function idw(z, dis, i, j, power::Float64)
    invdisp = 1.0 ./ (dis.^power)
    weights = invdisp ./ sum(invdisp)
    dot(weights, z)
end


function blocking_idw(a, nodata, radius, power::Float64)
    nrow, ncol = size(a)
    nnp = 12 # number of nearest points

    # previously used library, KDTrees, included points exactly on the radius
    # NearestNeighbors excludes those. To keep results the same, use prevfloat
    radius = prevfloat(Float64(radius))

    ai = copy(a) # copy for filling in
    # transform the array with missing data to XY and Z
    xy, z = to_xy_z(a, nodata)
    xyint = round.(Int, xy) # the integer version is for the line drawing performance
    tree = BallTree(xy) # construct a BallTree for fast inrange searches
    barriers = a .!= nodata

    p = Progress(ncol * nrow, 1, "Interpolating...")
    for j = 1:ncol
        for i = 1:nrow

            # skip data cells
            if a[i,j] != nodata
                next!(p)
                continue
            end

            # find indices of closest N points
            ind = inrange(tree, Float64[i, j], radius)

            # skip cells with no neighbours within radius
            if isempty(ind)
                ai[i,j] = nodata
                next!(p)
                continue
            end

            pivot = Float64[i, j]
            neighbours = xy[:, ind]
            dis = colwise(Euclidean(), pivot, neighbours)

            # find the barricaded neighbours
            idxs = Int[] # list of indices that are not barricaded (positive selection)
            sizehint!(idxs, nnp)
            for (idx, np) in enumerate(ind)
                if !isbarricaded(barriers, i, j, xyint[1, np], xyint[2, np])
                    push!(idxs, idx)
                end
            end

            # lose the barricaded neighbours
            ind_nonbarricaded = ind[idxs]
            dis_nonbarricaded = dis[idxs]

            ai[i, j] = idw(z[ind_nonbarricaded], dis_nonbarricaded, i, j, power)
            next!(p)
        end
    end
    ai
end

end # module
