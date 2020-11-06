#= 
filter Clouds based on attributes and surrounding points
Format should be according to:

    classify!cloud::Cloud; kwargs)
        cloud[subset]
    end

Dirk Eilander, Martijn Visser, Deltares, 11-2016 =#

"function returns cloud indices based on statistical reduce function "
function reduce_index(cloud::Cloud, r::Raster;
    reduceri=reducer_minz_index, # function to find index of min max
    pointfilter=nothing, # predicate function to filter individual points
    min_dens=0) # minimum point density in cell to consider for educed cloud)

    subset = Int32[]
    # loop through cells
    @showprogress 5 "Reducer processing.." for i in 1:length(r)
        # find points in range
        idx0 = r[i]
        # loop through points and filter
        if pointfilter != nothing # if filter is given
            idx = Int32[] #  indices after filter
            for i in idx0
                pointfilter(cloud, i) || continue
                push!(idx, i)
            end
        else
            idx = idx0 # no filter applied
        end

        np = length(idx)
        # calculate density [points / m2]
        if min_dens > 0
            di = np / area
            di < min_dens && continue # min density threshold
        end

        np == 0 && continue # stat functions don't work on empty arrays
        # set statistics to grid
        push!(subset, idx[reduceri(cloud, idx)])
    end
    subset
end

"function to reduce pointcloud based on statistical min per gridcell"
function reduce_min(cloud::Cloud, cellsize::Real;
    pointfilter=nothing, # predicate function to filter individual points
    min_dens=0) # minimum point density in cell to consider for reduced cloud

    r = define_raster(cloud, cellsize; pointfilter=pointfilter)
    reduce_min(cloud, r; pointfilter=pointfilter, min_dens=min_dens)
end

"function to reduce pointcloud based on min per gridcell"
function reduce_min(cloud::Cloud, r::Raster;
    pointfilter=nothing, # predicate function to filter individual points
    min_dens=0) # minimum point density in cell to consider for reduced cloud

    subset = reduce_index(cloud, r; reduceri=reducer_minz_index, pointfilter=pointfilter, min_dens=min_dens)

    # return raster with statistics, density is saved to last layer
    cloud[subset]
end

"function to reduce pointcloud based on statistical max per gridcell"
function reduce_max(cloud::Cloud, cellsize::Real;
    pointfilter=nothing, # predicate function to filter individual points
    min_dens=0) # minimum point density in cell to consider for educed cloud

    r = define_raster(cloud, cellsize; pointfilter=pointfilter)
    reduce_max(cloud, r; pointfilter=pointfilter, min_dens=min_dens)
end

"function to reduce pointcloud based on statistical max per gridcell"
function reduce_max(cloud::Cloud, r::Raster;
    pointfilter=nothing, # predicate function to filter individual points
    min_dens=0) # minimum point density in cell to consider for educed cloud

    subset = reduce_index(cloud, r; reduceri=reducer_maxz_index, pointfilter=pointfilter, min_dens=min_dens)

    # return raster with statistics, density is saved to last layer
    cloud[subset]
end

function reduce_pointfilter(cloud::Cloud, pointfilter)
    subset = Int32[]
    for i in 1:length(positions(cloud))
        (pointfilter != nothing) && (pointfilter(cloud, i) || continue)
        push!(subset, i)
    end
        
    cloud[subset]
end
