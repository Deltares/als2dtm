using GridOperations; const GO = GridOperations
using Base.Test
using ColorTypes

@testset "Colorize" begin
    arr = [3.0, 4.0, 5.0, NaN]
    colored = GO.colorize(arr)
    @test colored == [
        RGB{Float64}(0.267004, 0.004874, 0.329415),
        RGB{Float64}(0.1281485, 0.565107, 0.5508925),
        RGB{Float64}(0.993248, 0.906157, 0.143936),
        RGB{Float64}(1.0, 1.0, 1.0)
    ]
end

@testset "Progressive morphological filter" begin
    srand(1)
    zmin = rand(3.0f0:0.1f0:5.0f0, 30, 40)
    max_window_radius = 9.0
    slope = 0.151234
    dh_max = 2.0
    dh_min = 0.4
    nodata = -9999.0
    boundarymask = trues(zmin)
    zmax_pmf, flags = GO.pmf_filter(zmin, boundarymask, max_window_radius, slope, dh_max, dh_min, 1.0)
    @test count(flags) === 754
    @test minimum(zmax_pmf) ≈ 3.4
    @test maximum(zmax_pmf) ≈ 3.702468  # large Float32 inaccuracy, should be 3.7
    @test hash(zmax_pmf) === 0x520f02f7bbe47b0e
end
