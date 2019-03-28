#=
Predicate functions to filter points in Clouds based on positions or attributes.
A point is defined by a PointCould and index.

onsistent use of pointfilters:
value in scripts can either be <nothing> or a function according to format below
--> true: the points are used in the function.
--> False: the points are skipped

Format should be according to:

    function pointfilter(cloud::Cloud, index::Int)
        true
    end

=#

ground_water(cloud::Cloud, i::Integer) = (cloud[:classification][i] == 2) || (cloud[:classification][i] == 9)
ground(cloud::Cloud, i::Integer) = (cloud[:classification][i] == 2)
notground(cloud::Cloud, i::Integer) = (cloud[:classification][i] != 2)
water(cloud::Cloud, i::Integer) = (cloud[:classification][i] == 9)
notwater(cloud::Cloud, i::Integer) = (cloud[:classification][i] != 9)
notoutlier(cloud::Cloud, i::Integer) = cloud[:classification][i] != 7 && cloud[:classification][i] != 18
outlier(cloud::Cloud, i::Integer) = cloud[:classification][i] == 7 || cloud[:classification][i] == 18
unclassified(cloud::Cloud, i::Integer) = cloud[:classification][i] == 0
unclassified_lastreturn(cloud::Cloud, i::Integer) = lastreturn(cloud, i) && unclassified(cloud, i)
unclassified_firstreturn(cloud::Cloud, i::Integer) = firstreturn(cloud, i) && unclassified(cloud, i)

highscanangle(cloud::Cloud, i::Integer) = abs(cloud[:scan_angle][i]) > 20.0
lowscanangle(cloud::Cloud, i::Integer) = abs(cloud[:scan_angle][i]) <= 20.0

# dummy filters
filtertrue(::Cloud, ::Integer) = true
filterfalse(::Cloud, ::Integer) = false

# ASPRS, but only on normalized pointcloud(!)
low_veg(cloud::Cloud, i::Integer) = 0.5 < getz(cloud, i) <= 2.0
med_veg(cloud::Cloud, i::Integer) = 2.0 < getz(cloud, i) <= 5.0
high_veg(cloud::Cloud, i::Integer) =  5.0 < getz(cloud, i)

non_tree(cloud::Cloud, i::Integer) = getz(cloud, i) < 50.0  # higher trees don't exist ;)

# check if last return
lastreturn(cloud::Cloud, i::Integer) = cloud[:return_number][i] == cloud[:number_of_returns][i]
firstreturn(cloud::Cloud, i::Integer) = cloud[:return_number][i] == 1
singlereturn(cloud::Cloud, i::Integer) = cloud[:return_number][i] == cloud[:number_of_returns][i] == 1

# edge of flightline
flightedge(cloud::Cloud, i::Integer) = cloud[:edge_of_flight_line][i] == 1

# combined

"First returns, not ground and not on edge. Default for chm method"
function unclassifiedfirstreturnonedge(cloud::Cloud, i::Integer)
    unclassified_firstreturn(cloud, i) && !flightedge(cloud, i)
end

"First returns, not ground and not with high scan angle."
function unclassifiedfirstreturnlowscanangle(cloud::Cloud, i::Integer)
    unclassified_firstreturn(cloud, i) && lowscanangle(cloud, i)
end

"First returns, not ground, only trees and not with high scan angle."
function unclassifiedfirstreturnnontreelowscanangle(cloud::Cloud, i::Integer)
    unclassified_firstreturn(cloud, i) && lowscanangle(cloud, i) && non_tree(cloud, i)
end

"not outlier and is last return. default pointfilter for classify_ground! method"
function notoutlier_lastreturn(cloud::Cloud, i::Integer)
    lastreturn(cloud, i) && notoutlier(cloud, i)
end

"not outlier and is last return. default pointfilter for classify_ground! method"
function notoutlier_ground(cloud::Cloud, i::Integer)
    notoutlier(cloud, i) && ground(cloud, i)
end

"not outlier and is last return. default pointfilter for classify_ground! method"
function notoutlier_notwater_lastreturn(cloud::Cloud, i::Integer)
    lastreturn(cloud, i) && notoutlier(cloud, i) && notwater(cloud, i)
end
