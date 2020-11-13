using als2dtm

infile = "data/small.laz"
outfolder = "output/"
epsg = 32748

# PMF settings
radius = 16.
slope = 0.6
dhmax = 1.1
dhmin = 0.5

lidar_pipeline(infile, epsg, outfolder, radius, slope, dhmax, dhmin)
