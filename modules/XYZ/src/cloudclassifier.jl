#=
Classify Clouds based on attributes and surrounding points
Format should be according to:

    classify!cloud::Cloud; kwargs)
        cloud'
    end

=#

## BASICS
" classification of point in Cloud"
function classify!(cloud::Cloud, i::Integer, class::UInt8)
    cloud[:classification][i] = class
    nothing
end


" classification of of Cloud based on pointfilter"
function classify!(cloud::Cloud, pointfilter, class::UInt8)
    for i in 1:length(cloud)
        (pointfilter(cloud, i) || continue)
        cloud[:classification][i] = class
    end
    nothing
end

"set classification of a point to zero"
function reset_class!(cloud::Cloud, i::Integer)
    classify!(cloud, i, UInt8(0))
    nothing
end

"set classification of all points in Cloud to zero"
function reset_class!(cloud::Cloud;
    pointfilter = nothing)

    for i in 1:length(cloud)
      (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
      classify!(cloud, i, UInt8(0))
    end
    nothing
end

## Raster based
"function to classify minima in pointcloud based per gridcell"
function classify_min!(cloud::Cloud, cellsize::Real, class::Int;
   pointfilter = nothing, # predicate function to filter individual points
   min_dens = 0) # minimum point density in cell to consider for educed cloud

    r = define_raster(cloud, cellsize; pointfilter = pointfilter)
    classify_min!(cloud, r, class; pointfilter = pointfilter, min_dens = min_dens)
end

"function to classify minima in pointcloud based per gridcell"
function classify_min!(cloud::Cloud, r::Raster, class::Int;
    pointfilter = nothing, # predicate function to filter individual points
    min_dens = 0) # minimum point density in cell to consider for educed cloud

    subset = reduce_index(cloud, r; reduceri = reducer_minz_index, pointfilter = pointfilter, min_dens = min_dens)

    # return raster with statistics, density is saved to last layer
    for i in subset
        classify!(cloud, i, UInt8(class))
    end

    nothing
end

"function to classify maxima in pointcloud based per gridcell"
function classify_max!(cloud::Cloud, cellsize::Real, class::Int;
    pointfilter = nothing, # predicate function to filter individual points
    min_dens = 0) # minimum point density in cell to consider for educed cloud

    r = define_raster(cloud, cellsize; pointfilter = pointfilter)
    classify_max!(cloud, r, class; pointfilter = pointfilter, min_dens = min_dens)
end

"function to classify maxima in pointcloud based per gridcell"
function classify_max!(cloud::Cloud, r::Raster, class::Int;
    pointfilter = nothing, # predicate function to filter individual points
    min_dens = 0) # minimum point density in cell to consider for educed cloud

    subset = reduce_index(cloud, r; reduceri = reducer_maxz_index, pointfilter = pointfilter, min_dens = min_dens)

    # return raster with statistics, density is saved to last layer
    for i in subset
        classify!(cloud, i, UInt8(class))
    end

    nothing
end

"classifies points with a cell if z <= to the z-value in surface + tolerance"
function classify_below_surface!(
    cloud::Cloud, r::Raster, surf::AbstractArray{T, 2}, class::Integer;
    tolerance = 0.0,
    pointfilter = nothing,
    nodata = -9999.0) where T <: Real

    ((size(r)[1] == size(surf)[1]) && (size(r)[2] == size(surf)[2])) || error("raster and surface size mismatch")
    for icell in 1:length(r)
        z_max = surf[icell] + tolerance
        z_max == nodata && continue
        for i in r[icell] # indices Cloud
            positions(cloud)[i][3] <= z_max || continue
            (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
            classify!(cloud, i, UInt8(class))
        end
    end

    nothing
end

function classify_below_surface!(cloud::Cloud, surf::Real, class::Integer;
    tolerance = 0.0,
    pointfilter = nothing)

   z_max = surf + tolerance
   for i in 1:length(position(cloud)) # indices Cloud
      positions(cloud)[i][3] <= z_max || continue
      (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
      classify!(cloud, i, UInt8(class))
   end

   nothing
end

"classifies points with a cell if z > to the z-value in surface + buffer"
function classify_above_surface!(cloud::Cloud, r::Raster,
    surf::AbstractArray{T, 2}, class::Integer;
    buffer = 0.0,
    pointfilter = nothing,
    nodata = -9999.0) where T <: Real

    ((size(r)[1] == size(surf)[1]) && (size(r)[2] == size(surf)[2])) || error("raster and surface size mismatch")
    for icell in 1:length(r)
        z_min = surf[icell] + buffer
        z_min == nodata && continue
        for i in r[icell] # indices Cloud
          positions(cloud)[i][3] > z_min || continue
          (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
          classify!(cloud, i, UInt8(class))
        end
    end

    nothing
end

"classify points outside of raster mask where true indicates data not to classify"
function classify_mask!(cloud::Cloud, r::Raster, mask::Union{BitArray, Array{Bool}}, class::Integer)
    size(r) == size(mask) || error("raster and mask size mismatch")
    for icell in 1:length(r)
        mask[icell] && continue  # skip points inside mask
        for i in r[icell] # pointcloud indices
          classify!(cloud, i, UInt8(class))
        end
    end
    nothing
end

function classify_above_surface!(cloud::Cloud, surf::Real, class::Integer;
    buffer = 0.0,
    pointfilter = nothing)

    z_min = surf + buffer
    for i in 1:length(position(cloud)) # indices Cloud
      positions(cloud)[i][3] > z_min || continue
      (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
      classify!(cloud, i, UInt8(class))
    end

    nothing
end

## Slope based
"""Adapted from Vosselman (2000)
slope based classification of pointdata from point cloud"""
function vosselman(cloud::Cloud, i::Integer, neighbors::Vector{U};
    max_slope=0.2,
    tolerance=0.0) where U <: Integer

    if length(neighbors) > 1
        # get positions point i
        xi, yi, zi = positions(cloud)[i]

        # get z values
        z_vec = getz(cloud, neighbors) # z only
        idx_order = sortperm(z_vec) # loop through points from min to max elevation

        # initiate
        isclass = false
        for j in neighbors[idx_order]
            cloud[:classification][j] == 7 && continue # skip comparison with low outlier
            # 1) get horizontal and vertical distance
            xj, yj, zj = positions(cloud)[j]
            di = hypot(xi - xj, yi - yj)
            dz = zj - zi + tolerance
            # 2) calculate max z difference based on distance
            thresh = -1 * di * max_slope
            # 3) compare with actual z diff; if z diff too large -> no ground point
            dz < thresh && break
            # 4) classify if all points below point i have been checked without break
            if dz > 0
                isclass = true
                break
            end
        end
      else
          isclass = true # no neighbors
      end

    isclass
end

## applied
"classify outliers"
function classify_outliers!(cloud::Cloud;
    cellsize=100.0, # cellsize size for regural grid
    dz=1.0, # min distance for outlier classification
    max_outliers=15, # max number of outlier per cell
    max_height=Inf) # maximum heigh above lowest non outlier in cell

    low_noise = UInt8(7) # from ASPRS Standard LIDAR Point Classes
    high_noise = UInt8(18) # from ASPRS Standard LIDAR Point Classes (LAS 1.4)
    do_high_noise = !isinf(max_height)
    # all points below -5 m are always set to low noise
    const cutoff_low = -5.0
    too_low(cloud::Cloud, i::Integer) = getz(cloud, i) < cutoff_low
    classify!(cloud, too_low, low_noise)
    # all points above 100 m are always set to high noise
    const cutoff_high = 100.0
    too_high(cloud::Cloud, i::Integer) = getz(cloud, i) > cutoff_high
    classify!(cloud, too_high, high_noise)
    # points outside the cutoff should not be included in the outlier analysis
    not_cut_off(cloud::Cloud, i::Integer) = cutoff_low <= getz(cloud, i) <= cutoff_high
    r = define_raster(cloud, cellsize, pointfilter=not_cut_off)
    # create dz Vector with length of pointcloud to store as new attribute later
    dz_outliers = zeros(length(positions(cloud)))
    # loop through cells
    @showprogress 5 "Classifying outliers.." for icell in 1:length(r)
        # find points in range
        idxs = r[icell] # get points in rastercell
        np = length(idxs)
        np < (max_outliers*10) && continue # max 10% outliers. make sure <max_outliers> is in balance with <cellsize>

        # get z values of points in cell
        z_vec = getz(cloud, idxs)
        if np >= 1e5 # RadixSort only faster with many points
            idx_order = sortperm(z_vec, alg=RadixSort)
        else
            idx_order = sortperm(z_vec)
        end

        # get difference in z for lowest max_outliers+1 points
        dz_vec = z_vec[idx_order[2:max_outliers+1]] - z_vec[idx_order[1:max_outliers]]
        if sum(dz_vec .>= dz) == 0 && !do_high_noise
            continue # no low outliers
        end
        ilast = findlast(dz_vec .>= dz) # find the index of the last jump
        zmin = z_vec[idx_order[ilast+1]] # lowest non outlier

        if do_high_noise
            # classify points higher than max_height above lowest non outlier as high noise
            zdif = z_vec - zmin
            toohi = zdif .> max_height
            itoohi = idxs[toohi]
            dz_outliers[itoohi] = zdif[toohi]
            for i in itoohi
                classify!(cloud, i, high_noise)
            end
        end

        for i in idxs[idx_order[1:ilast]]
            dz_outliers[i] = zmin - getz(cloud, i) # difference in elevation with lowes non outlier
            classify!(cloud, i, low_noise)
        end
    end
    # update Cloud with dz of outliers
    cloud[:dz] = dz_outliers
    nothing
end

"Apply Vosselman (2000) to classify ground points"
function classify_ground!(cloud::Cloud, spatialindex::NearestNeighbors.NNTree;
    radius=3.0, # radius to check slope condition [m]
    max_slope=0.3, # max allowd slope between ground points [m/m]
    tolerance=0.0, # vertical tolerance [m]
    pointfilter=nothing, # only classify points with class <class_only>
    class=2)

    # check if the spatialindex is 2D Euclidean
    @assert (spatialindex.metric == Euclidean()) && (length(spatialindex.data[1]) == 2)
    n = length(positions(cloud))
    # loop through points
    @showprogress 5 "Classifying low surface points (Vosselman).." for i in 1:n
        # get xyz & class of point i
        (pointfilter != nothing) && (pointfilter(cloud, i) || continue)

        xi, yi, zi = positions(cloud)[i]
        # find number of points in range
        neighbors = inrange(spatialindex, [xi, yi], radius)

        # classify
        isground = vosselman(cloud, i, neighbors; max_slope=max_slope, tolerance=tolerance)

        # classes according to asprs las format
        if isground
            classify!(cloud, i, UInt8(class)) # ground
        end
    end

    nothing
end

"""Apply Vosselman (2000) to classify water points
only apply to points below min_z_tolerance from minima in grid with cellsize=radius"""
function classify_water!(cloud::Cloud, spatialindex::NearestNeighbors.NNTree;
    radius=50.0, # radius to check slope condition [m]
    max_slope=0.01, # max allowd slope between ground points [m/m]
    tolerance=0.05, # vertical tolerance vosselman [m] (due to vairance in water level measurements)
    min_z_tolerance = 0.3, # only apply to points below min_z_tolerance from minima in grid with cellsize=radius
    pointfilter = nothing, # only classify points with class <class_only>
    class = 9)

    r_large = define_raster(cloud, radius; pointfilter = pointfilter)
    z_max = rasterize(cloud, r_large, reducer = reducer_minz) + min_z_tolerance

    # check if the spatialindex is 2D Euclidean
    @assert (spatialindex.metric == Euclidean()) && (length(spatialindex.data[1]) == 2)

    # loop through points
    @showprogress 5 "Classifying water points (Vosselman).." for icell in 1:length(r_large)
        max_z = z_max[icell]
        idx0 = r_large[icell]
        for i in idx0
            # get xyz & class of point i
            (pointfilter != nothing) && (pointfilter(cloud, i) || continue)

            xi, yi, zi = positions(cloud)[i]
            zi <= max_z || continue

            # find number of points in range
            neighbors = inrange(spatialindex, [xi, yi], radius)

            # classify
            iswater = vosselman(cloud, Int(i), neighbors; max_slope=max_slope, tolerance=tolerance)

            # classes according to asprs las format
            if iswater
                classify!(cloud, i, UInt8(class))  # ground
            end
        end
    end

    nothing
end
