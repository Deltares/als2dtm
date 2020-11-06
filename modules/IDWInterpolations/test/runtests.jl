using IDWInterpolations
using Base.Test

function testdata(nrow, ncol, nodata)
    srand(0)
    a = rand(2.0:8.0, nrow, ncol)
    for i = 1:nrow * ncol
        if rand() > 0.2
            a[i] = nodata
        end
    end
    a
end

nodata = 0.0
a = testdata(4, 6, nodata)
ai = blocking_idw(a, nodata, 2, 2.0)

ai_check = [8.0 8.0 7.5 7.0 7.666666666666666 8;
    6.666666666666666 8.0 5.666666666666666 7.0 7.333333333333332 7.2;
    4.0 4.285714285714285 2.0 5.571428571428571 8.0 6.0;
    2.0 2.5 2.0 5.0 5.0 6.0]

@test ai ≈ ai_check
