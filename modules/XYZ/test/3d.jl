using XYZ
using StaticArrays
using BenchmarkTools, Compat

a = collect(1:16)
a = reshape(a, (2,2,4))
XYZ.offsets(a)

bbox = XYZ.NBoundingBox(xmin=542900.0, xmax=544100.0, ymin=9679900.0, ymax=9681077.6, zmin=1.06, zmax=329.09)
dims = SVector(5, 5, 5)
coords = MVector(1.0,2.0,3.0)
coords = MVector(1,-1,2)
@btime invalid_coords(dims, coords)
@btime sub2ind(dims, coords[1], coords[2], coords[3])
