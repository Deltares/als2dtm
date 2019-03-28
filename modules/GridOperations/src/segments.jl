#function for segmentation in grids based on thresholding and connectivity
#
# based on functions from Images module


function im_labeling(tf::Matrix{U}; connectivity=8) where U <: Integer
    # create connectivity array
    if connectivity == 8
        connect_array = trues(3,3)
    elseif connectivity == 4
        connect_array = BitMatrix([0 1 0; 1 1 1; 0 1 0])
    end

    # label areas
    segments = label_components(tf, connect_array)
end

# copy of function above, not sure how to catch both types
function im_labeling(tf::BitMatrix; connectivity=8)
    # create connectivity array
    if connectivity == 8
        connect_array = trues(3,3)
    elseif connectivity == 4
        connect_array = BitMatrix([0 1 0; 1 1 1; 0 1 0])
    end

    # label areas
    segments = label_components(tf, connect_array)
end

"""cleanup binary images based on minimal area of segments
Input
segments        labeled array
min_area        min area of one segment [m2]
res             resolution of 1 pixel [m2]
connectivity    all surrounding pixels (=8) or only N-E-S-W pixels (=4) (default = 8)

output
updated label array
"""
function im_segments_cleanup(segments::Array{U, 2}, min_area::Real;
        res=1, connectivity=8) where U <: Integer
    # remove noise
    segments = copy(segments)
    segment_areas = component_lengths(segments)
    segment_idx = component_indices(segments)
    for lab in 1:maximum(segments) # skip background == 0
        area = segment_areas[lab+1]
        if area < (min_area / (res*res))
            for idx in segment_idx[lab+1]
                segments[idx] = 0
            end
        end
    end
    # re-label segments
    im_labeling(segments; connectivity=connectivity)
end


function split_segments(segments::Array{U, 2}, block_size, min_area;
    res=1, connectivity=8) where U <: Integer

    # split segments
    segments = copy(segments)
    segment_idx = component_subscripts(segments)
    segment_bbox = component_boxes(segments)
    block_size = Float64(block_size) / res
    nsegments = maximum(segments)
    for lab in 1:nsegments # skip background == 0
        (ymin, xmin), (ymax, xmax) = segment_bbox[lab+1]
        dx, dy = xmax-xmin+1, ymax-ymin+1
        if ((dx > 1.5*block_size) || (dy > 1.5*block_size))
            nrow, ncol = min(Int(cld(dy, block_size)),1), min(Int(cld(dx,block_size)),1)
            for (y, x) in segment_idx[lab+1]
                # cells include bottom and left boundaries
                col = Int(fld(x - xmin, block_size)+1)
                row = Int(fld(ymax - y, block_size)+1)
                iblock = sub2ind((nrow, ncol), row, col)
                segments[y,x] = iblock
            end
        end
    end
    # re-label segments
    segments = im_labeling(segments; connectivity=connectivity)

    # merge segments with too small area
    segment_areas = component_lengths(segments)
    segment_idx = component_indices(segments)
    for lab in 1:maximum(segments) # skip background == 0
        area = segment_areas[lab+1]
        if area < (min_area / (res*res))
            b = falses(segments)
            for idx in segment_idx[lab+1]
                b[idx] = true
            end
            bd = dilate(b)
            neighboring_cells = bd .& broadcast(~, b)
            new_lab = maximum(segments[neighboring_cells]) # find label neighboring segment
            new_lab == 0 && continue # no neighboring segments
            for idx in segment_idx[lab+1]
                segments[idx] = new_lab
            end
        end
    end
    # re-label segments
    im_labeling(segments; connectivity=connectivity)
end

"""reduce pixel in 'A' within a connected segment in binary grid 'tf' using function 'func'

Input
segments        labeled array
A               input image with to sample from can [Array{2,Real}]

output
image with reduced value in all segment pixels"""
function im_segments_reduce(segments::Array{U, 2}, A::Array{T, 2};
                            connectivity=8,
                            func=x -> maximum(x)-minimum(x),
                            nodata=-9999.0) where U <: Integer where T <: AbstractFloat
    @assert size(segments) == size(A) "the size of image A and binary grid should match"

    # get indices
    initiate = true
    B = 0

    for idx in component_indices(segments)
        if initiate # skip background, but initiate output
            # initiate output array with output type at first loop
            vals = func(T[1, 2])
            @assert isa(vals, Real) "incorrect reducer function"
            reduced_type = typeof(vals)
            B = fill(reduced_type(nodata), size(A))
            initiate = false
        else
            vals = [xi for xi=A[idx] if xi != nodata]
            if !isempty(vals)
                B[idx] = func(vals)
            end
        end
    end

    B
end
