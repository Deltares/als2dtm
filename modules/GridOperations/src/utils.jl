
function ind2sub(a, i)
    i2s = CartesianIndices(a)
    i2s[i].I
end
function sub2ind(a, i...)
    s2i = LinearIndices(a)
    s2i[i...]
end

"""tests if point p1 and p2 are on the same side of line ab"""
function same_side(p1::Vector{U}, p2::Vector{U}, a::Vector{U}, b::Vector{U}) where U <: AbstractFloat
    cp1 = cross(b - a, p1 - a)
    cp2 = cross(b - a, p2 - a)
    dot(cp1, cp2) >= 0
end

""" tests if a point p is Left|On|Right of an infinite line ab
- >0 for p left of the line ab
- =0 for p on the line
- <0 for p right of the line
"""
function is_left(p::Vector{U}, a::Vector{U}, b::Vector{U}) where U <: Real
    (b[1] - a[1]) * (p[2] - a[2]) - (p[1] - a[1]) * (b[2] - a[2])
end

"""calculates 2D convec hull for set of points
returns points which determine the hull and its indices in the input list

based on Andrew's Monotone Chain Algorithm (Andrew, 1979), explained here:
http://geomalgorithms.com/a10-_hull-1.html"""
function get_hull(P::Vector{Vector{T}}) where T <: Real
    np = length(P)
    # sort indices, first by x then y
    pm = zeros(T, (np, 3))
    [pm[i,:] = [P[i][1], P[i][2], i] for i in 1:np] # to 2d array
    pm = sortrows(pm, by=x -> (x[1], x[2]))
    P = [pm[i,1:2] for i in 1:np]
    idx = Int[pm[i,3] for i in 1:np]
    # get indices of points with 1st x min or max and 2nd y min or max
    min_min = 1 # min x, min y
    i = 0 ## initiate so it is known outside loop, arghh...
    for i in 2:np
        P[i][1] != P[1][1] && break
    end
    min_max = i - 1
    for i in np - 1:-1:1
        P[i][1] != P[np][1] && break
    end
    max_min = i + 1
    max_max = np
    # initiate convex hull stack with point indices
    stack = Int[]

    ## LOWER HULL
    push!(stack, min_min)
    top = 1 # no. of points in stack
    for i in 2:max_min - 1 # loop through points with increasing x sequence
        is_left(P[i], P[min_min], P[max_min]) >= 0 && continue # ignore P[i] above or on the lower line
        while top >= 2 # at least two points in stack
            is_left(P[i], P[stack[top - 1]], P[stack[top]]) > 0 && break # P[i] in hull
            pop!(stack)
            top -= 1
        end
        push!(stack, i)
        top += 1
    end
    push!(stack, max_min)
    top += 1

    ## UPPER HULL
    # if  distinct xmax points push max_min point onto stack
    if max_max != max_min
        push!(stack, max_max)
        top += 1
    end
    bot = top # bottom point of upper hull
    for i in np - 1:-1:2 # loop through points with decreasing x sequence
        is_left(P[i], P[max_max], P[min_max]) >= 0 && continue # ignore P[i] below or on the upper line
        while bot > top # at least two points in upper stack
            is_left(P[i], P[stack[top - 1]], P[stack[top]]) > 0 && break # P[i] in hull
            pop!(stack)
            top -= 1
        end
        push!(stack, i)
        top += 1
    end
    push!(stack, min_max)
    top += 1

    # return points of hull and indices
    P[stack], idx[stack]
end

"""calculates bbox coordinates for set of points
returns the list with coordinates the indices of the points which determine the bbox"""
function get_bbox(P::Vector{Vector{T}}) where T <: Real
    @assert length(P) >= 2
    xmin, xmax, ymin, ymax = Inf, -Inf, Inf, -Inf
    xmini, xmaxi, ymini, ymaxi = 0, 0, 0, 0  # declare indices of points that make up bbox
    i = 0
    @inbounds for p in P
        addp = false
        i += i  # index of point in vector P
        if p[1] < xmin
            xmin = p[1]
            xmini = i
        end
        if p[1] > xmax
            xmax = p[1]
            xmaxi = i
        end
        if p[2] < ymin
            ymin = p[2]
            ymini = i
        end
        if p[2] > ymax
            ymax = p[2]
            ymaxi = p
        end
    end

    # return bbox and indices that make up bbox
    [xmin,xmax,ymin,ymax], [xmini, xmaxi, ymini, ymaxi]
end

function triangle_area(a::Vector{U}, b::Vector{U}, c::Vector{U}) where U <: AbstractFloat
    abs(a[1] * (b[2] - c[2]) + b[1] * (c[2] - a[2]) + c[1] * (a[2] - b[2])) / 2.
end

"""check if point p [x,y] in triangle made out of points a,b,c"""
function in_triangle(p::Vector{U}, a::Vector{U}, b::Vector{U}, c::Vector{U}) where U <: AbstractFloat
    p, a, b, c = [[x[1],x[2],0.] for x in [p,a,b,c]]
    same_side(p, a, b, c) && same_side(p, b, a, c) && same_side(p, c, a, b)
end

"""check if point (x,y) in traingle made out of points vector"""
function in_triangle(x::AbstractFloat, y::AbstractFloat, points::Vector{Vector{U}}) where U <: AbstractFloat
    length(points)
    npoints != 3 && @warn("more than three points in list, check based on first three points only")
    a, b, c = points[1], points[2], points[3]
    xypoint_in_triangle([x, y], a, b, c)
end

"""check if point (x,y) in bbox"""
function in_bbox(x::AbstractFloat, y::AbstractFloat, points::Vector{Vector{U}}) where U <: AbstractFloat
    xmin, xmax, ymin, ymax = get_bbox(points)[1]
    xmin <= x <= xmax && ymin <= y <= ymax
end


"""check if point in hull of all points based on hull computation"""
function in_hull(x::Real, y::Real, points::Vector{Vector{U}}) where U <: Real
    npoints = length(points)
    @assert npoints >= 3 ["point vector should contain at least three points"]
    push!(points, [x,y])
    p_hull, idx_hull = get_hull(points)
    !in(npoints + 1, idx_hull)
end
