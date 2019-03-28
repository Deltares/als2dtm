"""percentile filter based on updating histograms while sliding over the image
implementation of Perreault and Hebert (2007) algorithm see: http://nomis80.org/ctmf.pdf

precision is limited by number of bins.
for a images with a bit depth of 8 bit (e.g. greyscale images), full precision
can be acchieved with a dual-level histogram of 16 bins
For images with a higher bit depth more bins are needed and the algorithm slows,
the precision can be set with the optional precision paramter

input
A           2d image Array
p           percentile value [1-100]
r           radius of window [n. cells]: i.e.: the window size is 2r+1 by 2r+1
is8bit      true if A has a bit-depth of 8 bits and no scaling needs to be applied (default false)
precision   number of decimals to preserve when binning
"""
function hist_filter(A::Array{T,2}, perc::Real, r::Int;
                                precision=0,
                                nodata_value=T(-9999.0),
                                hist_min=nothing,
                                hist_max=nothing) where T <: Real
    scale = T(10^precision) # default scale value
    ksize = (2*r+1)^2 # kernel size
    nrows, ncols = size(A) # image size
    mask = (A .≠ nodata_value) .& isfinite.(A)
    if T == UInt8
        im_min, im_max = 0, 255
    elseif T <: Integer
        precision > 0 && info("no scaling applied as A has integer values")
        scale = T(1)
        im_min, im_max = floor(Int, minimum(A[mask])), ceil(Int, maximum(A[mask]))
    else
        im_min, im_max = floor(Int, minimum(A[mask]))*scale, ceil(Int, maximum(A[mask]))*scale
    end
    # scale min max of hist to hist_min or hist_max if given
    if hist_min != nothing
        im_min = hist_min * scale
    end
    if hist_max != nothing
        im_max = hist_max * scale
    end

    nodata_value_scaled = nodata_value*scale
    @assert (im_max - im_min) >= 1 "the image bit depth with given 'precision', 'nodata_value', 'hist_min' and 'hist_max' params is zero"
    # create histogram bin size
    nbins = ceil(Int, sqrt(im_max-im_min))
    # fine level histograms have bins with size 1
    # course bin histogram bins are defined by:
    edges_course = Int[im_min-1:nbins:nbins*nbins+im_min;]

    # initialization of histograms
    # the counts arrays are not setup as offsetarrays because slicing doesn't work on these arrays
    # instead the arrays have an offset of r in indices
    cnts_course = zeros(Int32, (nrows+2*r, nbins))
    cnts_fine = zeros(Int32, (nrows+2*r, nbins, nbins))
    # create seperate bins for both extremes and nodata
    cnts_nodata = zeros(Int32, nrows+2*r)
    cnts_min = zeros(Int32, nrows+2*r)
    cnts_max = zeros(Int32, nrows+2*r)

    # create output array B
    B = similar(A)
    # pad A
    Apad = padarray(A, Pad(:reflect, r, r))

    function add_value!(i::Int, v::Real; w=1)
        if scale > 1
            v = round(v * scale, 0) # scaling for non 8bit images
        end
        if v ≈ nodata_value_scaled
            cnts_nodata[i+r] += w # offset r in cnts arrays
        elseif v < im_min
            cnts_min[i+r] += w
        elseif v >= im_max
            cnts_max[i+r] += w
        else
            icourse = min(Int(div(v-im_min, nbins)) + 1, nbins) # min to include right boundary
            ifine = min(Int(div(v-edges_course[icourse], 1)), nbins) # fine bins have size 1
            cnts_course[i+r, icourse] += w # offset r in cnts arrays
            cnts_fine[i+r, ifine, icourse] += w
        end

        nothing
    end

    # find bin corresponding to percentile start from FIRST bin
    function bin_count_course(i::Int, p_idx::Int, pcum::Int)
        pbin, icourse = 0, 0
        for icourse in 1:nbins
            # use of sum is faster than loop. tested with BenchmarkTools and ProfileView
            pbin = sum(cnts_course[i:i+r*2, icourse])
            pcum += pbin #[icourse]
            pcum >= p_idx && break
        end
        pcum -= pbin
        pcum, icourse
    end
    function bin_count_fine(i::Int, p_idx::Int, pcum::Int, icourse::Int)
        ifine=0
        for ifine in 1:nbins
            pcum += sum(cnts_fine[i:i+r*2, ifine, icourse])
            pcum >= p_idx && break
        end
        pcum, ifine
    end

    # find bin corresponding to percentile -> start from LAST bin
    function bin_count_course_rev(i::Int, p_idx::Int, pcum::Int)
        pbin, icourse = 0, 0
        for icourse in nbins:-1:1
            pbin = sum(cnts_course[i:i+r*2, icourse])
            pcum += pbin #[icourse]
            pcum >= (ksize - p_idx) && break
        end
        pcum -= pbin
        pcum, icourse
    end
    function bin_count_fine_rev(i::Int, p_idx::Int, pcum::Int, icourse::Int)
        ifine = 0
        for ifine in nbins:-1:1
            pcum += sum(cnts_fine[i:i+r*2, ifine, icourse])
            pcum >= (ksize - p_idx) && break
        end
        pcum, ifine
    end

    # find percentile value from multi level histogram
    function calc_perc(i::Int)
        # number of nodata_size values in kernel
        nodata_size = sum(cnts_nodata[i:i+r*2])
        if nodata_size < ksize # calculate percentile if more than half window with values
            p_idx = ceil(Int, max((ksize-nodata_size) * perc / 100., 1))
            # find percentile bin
            icourse, ifine = 0, 0
            if perc <= 100
                hist_min != nothing ? pcum = sum(cnts_min[i:i+r*2]) : pcum = 0
                # start from begin of histogram
                pcum, icourse = bin_count_course(i, p_idx, pcum)
                pcum, ifine = bin_count_fine(i, p_idx, pcum, icourse)
            else
                hist_max != nothing ? pcum = sum(cnts_max[i:i+r*2]) : pcum = 0
                # start from end in histograms
                pcum, icourse = bin_count_course_rev(i, p_idx, pcum)
                pcum, ifine = bin_count_fine_rev(i, p_idx, pcum, icourse)
            end

            # calculate the value of fine bin at
            return (edges_course[icourse] + ifine) / scale
        else
            return nodata_value
        end
    end

    ## start actual looping
    for j = 1:ncols, i in 1:nrows
        if j == 1 && i == 1 # initiate top left
            for kj = -r:r, ki = -r:r
                add_value!(i+ki, Apad[i+ki, j+kj])
            end
        elseif j == 1 && i > 1 # initiate first col
            for kj in -r:r
                add_value!(i+r, Apad[i+r, j+kj])
            end
        elseif i == 1 && j > 1 # first row
            for ki in -r:r
                add_value!(i+ki, Apad[i+ki, j-r-1]; w=-1) # remove old value
                add_value!(i+ki, Apad[i+ki, j+r])
            end
        else
            add_value!(i+r, Apad[i+r, j-r-1]; w=-1) # remove old value
            add_value!(i+r, Apad[i+r, j+r])
        end
        B[i,j] = calc_perc(i)
    end
    B
end


"""apply vosselman criterium to grid,

vosselman criterium:
z_{i} <=  min z[ki] + dist[ki] * max_slope
for ki in i-r:i+r, j-r:j+r

output
flags   true if criterium not met
B       max allow value per cell given neighbors
"""
function vosselman_filter(A::Array{T,2},
    radius,             # kernel radius [m]
    max_slope,          # mas tolorated slope [m/m]
    tolerance=0;
    res=1.,             # resolution of grid a.k.a. cell size [m]
    nodata::T=T(-9999.0)) where T <: AbstractFloat

    r = round(Int, radius/res) # length in # off cells in which to search for neighbours

    # create output array B
    nrows, ncols = size(A)
    B = ones(T, (nrows, ncols)) * Inf16
    flags = falses(A)
    # pad A
    Apad = ImageFiltering.padarray(A, Pad(:replicate, r, r))

    # loop offsets
    for kj = -r:r, ki = -r:r
        k_dist = hypot(ki,kj)
        ((k_dist > r) || ((ki==1) && (kj==1))) && continue # skip offset if out of range or zero
        k_tol = k_dist .* res .* max_slope + tolerance
        ## loop through image
        for j = 1:ncols, i = 1:nrows
            zi = A[i,j]
            zi == nodata && continue
            # replace if smaller than current (find min in kernel)
            zj = Apad[i+ki,j+kj] + k_tol
            zj == nodata && continue  # skip nodata values in kernel
            if zj < B[i, j]
                B[i, j] = zj
            end
            if zj < zi
                flags[i, j] = true
            end
        end
    end

    B, flags
end

# morphological opening & tophat = img - opening(img)
# square window, no circular kernel yet available!
function morph_tophat(img::Matrix{Float32}, boundarymask::BitMatrix, r::Int)
    # img_min = mapwindow(minimum, img, (s, s)) # erosion
    # img_open = mapwindow(maximum, img_min, (s, s)) # dilation
    img_min = morph_erosion(img, boundarymask, r)
    img_open = morph_dilation(img_min, boundarymask, r)
    img_open
end

"""Implementation of the Progressive Morphological Filter (Zhang, 2003)
works based on zmin grid where dropouts are filled with nearest neighbor value
NOTE: make sure A has no nodata values

output
flags   true if criterium not met
B       max allow value per cell given neighbors

(Zhang, 2003): http://users.cis.fiu.edu/~chens/PDF/TGRS.pdf """
function pmf_filter(A::Matrix{Float32}, boundarymask::BitMatrix,
    max_window_radius::Float64, # max radius for squared window [m]
    slope::Float64,             # terrain slope [m/m]
    dh_max::Float64,            # maximum elevation threshold between subsequent tophat operations [m]
    dh_min::Float64,            # initial elevation threshold [m]
    res::Float64)               # cellsize [m]

    @show max_window_radius, slope, dh_min, dh_max, res

    nrow, ncol = size(A)
    totcells = nrow * ncol
    # maximum window size of kernel in pixels
    rad_pix_max = round(Int, max_window_radius / res)
    # exponential increase radius; window size = 3,5,9,17,33,..,N
    kmax = ceil(Int, log2(rad_pix_max*2))
    w = [Int(exp2(k)) + 1 for k in 1:kmax]
    w[end] = rad_pix_max * 2 + 1
    # deprecate linear window growth -> does not give satisfying results
    # linear increase radius; window size = 3,5,7,9...,N
    # w = Int[2*k + 1 for k in 1:rad_pix_max]

    # give a warning when the slope will exceed dh_max
    slope * (w[end] - w[end-1]) * res + dh_min > dh_max && warn("dh_max will limit pmf filter")

    # initialize
    dh_t = dh_min
    Af = copy(A)
    # create max allowed surface (use to classify pointcloud)
    B = Af + dh_t
    # loop through kernel sizes

    @showprogress 1 "Progessive filtering..." for (k, wk) in enumerate(w)

        # threshold
        if k > 1
            dh_t = min(slope * (wk - w[k-1]) * res + dh_min, dh_max)
        end

        # calculate morph opening
        Af = morph_tophat(Af, boundarymask, fld(wk, 2))

        # update max allowed surface
        if k == 1
            B[:,:] = Af + dh_t
        else
            for I in eachindex(Af)
                # replace if smaller than current z_max
                z_max = Af[I] + dh_t
                B[I] = min(B[I], z_max)
            end
        end

    end
    # flag non ground points based on dh_t threshold
    flags = Array{Bool}(A .> B)
    # ncells = count(x -> x == 1, flags)
    # debugging print statement to be placed in loop above
    # println("  iter $(k): $(ncells)/$(totcells) filtered, window size = $(wk) cells, dh = $(@sprintf("%.2f", dh_t))")

    B, flags
end

# lee filter
# http://stackoverflow.com/questions/39785970/speckle-lee-filter-in-python
# Lee, J. S. (1980). Digital image enhancement and noise filtering by use of
# local statistics. IEEE transactions on pattern analysis and machine intelligence, (2), 165-168.
#
# enl (equivalen number of looks) -> 1 in LiDAR intensity data
# addititive noise
function lee_filter(img::Array{T,2}, size::Real; nodata=T(-9999.0)) where T <: AbstractFloat
    # local statistics
    fmean(A::AbstractArray) = mean(filter(x -> x != nodata, A))
    window_mean = mapwindow(fmean, img, (size, size))
    window_sqr_mean = mapwindow(fmean, img.^2, (size, size))
    window_variance = window_sqr_mean - window_mean.^2
    # global variance
    img_variance = var(filter(x -> x != nodata, img))


    # smoothing weights
    w = window_variance.^2 ./ (window_variance.^2 + img_variance.^2)

    # smoothed image
    img_output = (img .* w) + (window_mean .* (1 - w))
end

# enhanced lee filter
#
# Lopes, A., Touzi, R., & Nezry, E. (1990). Adaptive speckle filters and scene
# heterogeneity. IEEE transactions on Geoscience and Remote Sensing, 28(6), 992-1000.
# enl (equivalent number of looks) -> 1.0 for lidar data
function enhanced_lee_filter(img::Array{T,2}, size::Real;
        enl=1.0,
        k=1.0,
        nodata=T(-9999.0)) where T <: AbstractFloat

    max_cof_var = sqrt(1+2/enl)
    cu = 0.523/sqrt(enl)

    # local statistics
    fmean(A::AbstractArray) = mean(filter(x -> x != nodata, A))
    window_mean = mapwindow(fmean, img, (size, size))
    window_sqr_mean = mapwindow(fmean, img.^2, (size, size))
    window_variance = window_sqr_mean - window_mean.^2
    # in case of img being Float32 and having very high values,
    # such as sometimes happens with intensity
    # slightly negative values for window_variance can arise
    # these are here set to 0.0
    window_variance = clamp.(window_variance, 0.0, Inf)
    window_cof_var = sqrt.(window_variance) ./ window_mean
    window_cof_var[window_cof_var .== 0] = 0.001 # avoid zero division

    # smoothing weights
    # filter in heterogeneous area
    w = T.(exp.(-k .* (window_cof_var - cu) ./ (max_cof_var - window_cof_var)))
    # preserver original value for "point target" (e.g around edges)
    w[window_cof_var .>= max_cof_var] = 0.0
    # return local mean in homogeneous areas
    w[window_cof_var .<= cu] = 1.0

    # smoothed image
    img_output = (img .* (1-w)) + (window_mean .* w)
end

# morphological operations with circular kernel
function morph_erosion(A::Matrix{Float32}, boundarymask::BitMatrix, r::Int)
    nrows, ncols = size(A) # image size
    B = copy(A) # create output array B
    Apad = padarray(A, Pad(:replicate, r, r)) # pad A

    # loop offsets
    @inbounds for kj = -r:r, ki = -r:r
        ((hypot(ki,kj) > r) || ((ki==0) && (kj==0))) && continue # skip offset if out of range or zero
        ## loop through image
        for j = 1:ncols, i = 1:nrows
            boundarymask[i,j] || continue
            # replace if smaller than current (find min in kernel)
            zo = Apad[i+ki,j+kj]
            if zo < B[i, j]
                B[i, j] = zo
            end
        end
    end
    B
end

function morph_dilation(A::Matrix{Float32}, boundarymask::BitMatrix, r::Int)
    nrows, ncols = size(A) # image size
    B = copy(A) # create output array B
    Apad = padarray(A, Pad(:replicate, r, r)) # pad A

    # loop offsets
    @inbounds for kj = -r:r, ki = -r:r
        ((hypot(ki,kj) > r) || ((ki==0) && (kj==0))) && continue # skip offset if out of range or zero
        ## loop through image
        for j = 1:ncols, i = 1:nrows
            boundarymask[i,j] || continue
            # replace if larger than current (find max in kernel)
            zo = Apad[i+ki,j+kj]
            if zo > B[i, j]
                B[i, j] = zo
            end
        end
    end
    B
end
