"""interpolate missing values (nodata) in a grid
methods available:
nearest neighbor
ordinary kriging
inverste distance interpolation

inputs
img          the image to be interpolated
nodata       nodata value in img
mask         boolean array. only interpolate where true
mode         "nn", "kriging", "idw"
k            number of neighbors (default=10; for kriging and idw only)
p            idw power (default=10.0; for idw only)
extrapolate  check if point inside triangle of neighbors (default=false; for kriging and idw only)
"""
function interp_missing!(
    img::Array{T,2},
    nodata::T;
    mask=trues(size(img)),
    mode="nn",
    k=10, # number of neighboring points used in kringing and idw
    p=2.0, # idw power
    extrapolate=true) where T <: AbstractFloat

    # initiate interpolators or spatial index
    xy_p, z_p = xy_z_grid(img, nodata=nodata)

    # skip images with only no-data values
    if length(xy_p) > 0
        if mode in ["nn", "kriging", "idw"]
            tree = KDTree(xy_p)
        else
            error("mode $(mode) not implemented")
        end

        # interpolation
        ismissing = (img .== nodata) .& mask
        dims = size(img)
        for I in 1:length(img)
            ismissing[I] || continue
            x,y = ind2sub(dims, I)
            if mode == "nn"
                img[I] = interpolator_nn(x, y, xy_p, z_p, tree)
            elseif mode == "kriging"
                img[I] = interpolator_kriging(x, y, xy_p, z_p, tree; k=k, extrapolate=extrapolate)
            elseif mode == "idw"
                img[I] = interpolator_idw(x, y, xy_p, z_p, tree; p=p, k=k, extrapolate=extrapolate)
            end
        end
    end

    img
end

"""interpolate z at xy based on values z_p at xy_p"""
function interp2d(xy::Matrix{U}, xy_p::Matrix{U}, z_p::Vector{T};
    mode="nn",
    k=10, # number of neighboring points used in kringing and idw
    p=2.0, # idw power
    extrapolate=true,
    tree::NearestNeighbors.NNTree=KDTree(xy_p)) where T <: AbstractFloat where U <: Real

    # output
    n = size(xy,2)
    z = zeros(T,length(xy))

    # interpolation
    for I in 1:n
        x,y = xy[:,I]
        if mode == "nn"
            z[I] = interpolator_nn(x, y, xy_p, z_p, tree)
        elseif mode == "kriging"
            z[I] = interpolator_kriging(x, y, xy_p, z_p, tree; k=k, extrapolate=extrapolate)
        elseif mode == "idw"
            z[I] = interpolator_idw(x, y, xy_p, z_p, tree; p=p, k=k, extrapolate=extrapolate)
        end
    end

    z
end

function interpolator_kriging(x::Real, y::Real,
    xy_p::Array{U,2}, z_p::Array{T,1},
    tree::NearestNeighbors.NNTree=KDTree(xy_p);
    k=10, # number of neighboring points used in kriging and idw
    extrapolate=true,
    nodata=-9999.0) where T <: AbstractFloat where U <: Real

    idxs, dists = knn(tree, [x, y], k, true)
    if !extrapolate && !in_hull(x, y, [xy_p[:,i] for i in idxs])
        T(nodata)
    elseif length(idxs) > 0
        # define a covariance model
        γ = GaussianVariogram(sill=1.0f0, range=1.0f0, nugget=0.0f0)
        # define an estimator (i.e. build the Kriging system)
        ordkrig = OrdinaryKriging(xy_p[:,idxs], z_p[idxs], γ)
        # estimate at target location
        μ, _ = estimate(ordkrig, T[x,y]) # not using σ²
        μ
    else
        error("no points to interpolate")
    end
end

function interpolator_idw(x::Real, y::Real,
    xy_p::Array{U,2}, z_p::Array{T,1},
    tree::NearestNeighbors.NNTree=KDTree(xy_p);
    p=2.0, # idw power
    k=10, # number of neighboring points used in kringing and idw
    extrapolate=true,
    nodata=-9999.0) where T <: AbstractFloat where U <: Real

    idxs, dists = knn(tree, [x, y], k, true)
    if !extrapolate && !in_hull(x, y, [xy_p[:,i] for i in idxs])
        T(nodata)
    elseif length(idxs) > 0
        idw_interpolation(z_p[idxs], dists; power=p)
    else
        error("no points to interpolate")
    end
end

function interpolator_nn(x::Real, y::Real,
    xy_p::Array{U,2}, z_p::Array{T,1},
    tree::NearestNeighbors.NNTree=KDTree(xy_p)) where T <: AbstractFloat where U <: Real

    idx, dists = knn(tree, Float32[x, y], 1)
    z_p[idx[1]]
end


"""loop over image to find missings & make arrays with datapoints"""
function xy_z_grid(img::Array{T,2}; nodata::Real=-9999.0) where T
    n = count(x -> x != nodata, img)
    xy = zeros(Float32, (2, n))
    z = zeros(Float32, n)
    dims = size(img)
    i = 0
    for I in 1:length(img)
        if img[I] != nodata
            i += 1
            x,y = ind2sub(dims, I)
            # build xy and z arrays for KDTree
            xy[1, i] = x
            xy[2, i] = y
            z[i] = img[I]
        end
    end
    xy, z
end


"fill las dropouts, but keep no data areas next to strips open"
function create_mask(img::Array; nodata::Real=-9999.0, n=10)
    mask = img .!= nodata
    for i in 1:n
        ImageMorphology.dilate!(mask)
    end
    for i in 1:n
        ImageMorphology.erode!(mask)
    end
    Array{Bool,2}(mask)
end

"""calculate the Inverse Distance Weighting interpolation for one point"""
function idw_interpolation(z::Vector{U}, dis::Vector{U}; power=10.) where U <: AbstractFloat
    invdisp = 1.0 ./ (dis .^ power)
    weights = invdisp ./ sum(invdisp)
    dot(weights, z)
end

"""linear interpolation based fox xy withing triangle in 3D plane"""
function linear_interpolation(x::Real, y::Real, a::Vector{U},
        b::Vector{U}, c::Vector{U}) where U <: AbstractFloat

    t1 = ((b[1]-a[1])*(c[3]-a[3])-(c[1]-a[1])*(b[3]-a[3]))
    t2 = ((b[2]-a[2])*(c[3]-a[3])-(c[2]-a[2])*(b[3]-a[3]))
    n = ((b[1]-a[1])*(c[2]-a[2])-(c[1]-a[1])*(b[2]-a[2]))
    z = a[3] + t1/n*(y-a[2]) - t2/n*(x-a[1])
end
