# functions to colorize numeric arrays or grayscale images with color maps

"""Colorize an array by applying a colormap from `ColorSchemes`

For available colormaps see [this overview](http://juliagraphics.github.io/ColorSchemes.jl/stable/assets/figures/colorschemes.png).
"""
function colorize(A::AbstractArray;
    contrast::Bool=false, # enhance contrast by cutting the 2% on each side, like the QGIS default
    cmap::Symbol=:viridis, # matplotlib default
    nodata=-9999.0,
    vmin=nothing, # costly to filter twice
    vmax=nothing,
    bg=colorant"black") # black, color of nodata

    novalues = (A .== nodata) .| isnan.(A)
    values = broadcast(~, novalues)

    # calculate vmin and vmax
    # taking into account completely empty arrays
    if vmin == nothing || vmax == nothing
        if all(novalues)
            # return black image
            # shortcut because cannot reduce over empty collection
            return fill(bg, size(A))
        else
            # do it here to prevent passing twice over the data
            # if both vmin and vmax are nothing
            vmin_calculated, vmax_calculated = if contrast
                quantile(A[values], [0.02, 0.98])
            else
                extrema(A[values])
            end
        end
    end
    if vmin == nothing
        vmin = vmin_calculated
    end
    if vmax == nothing
        vmax = vmax_calculated
    end

    # get the right colormap from ColorSchemes
    colormap = getfield(ColorSchemes, cmap)
    if (vmin ≈ vmax) || (vmin > vmax)
        return fill(bg, size(A))
    end
    sc_img = (A .- vmin) .* (1.0 / (vmax - vmin))
    sc_img[novalues] = 0.0  # set to background color later, just to avoid get NaN InexactError
    color(x) = get(colormap, x)
    cimg = color.(sc_img)
    cimg[novalues] = bg
    cimg
end

colorize(A::Union{BitMatrix,Matrix{Bool}}) = Gray.(A)
