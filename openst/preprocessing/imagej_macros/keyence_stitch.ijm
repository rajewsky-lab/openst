arglist = getArgument();
args = split(arglist,';');
var gridx = args[0];
var gridy = args[1];
var indir = args[2];
var outfile = args[3];

run("Grid/Collection stitching", "type=[Grid: snake by rows] order=[Right & Down] grid_size_x=" + gridx + " grid_size_y=" + gridy + " tile_overlap=20 first_file_index_i=1 directory=[" + indir + "] file_names=Image_{iiiii}_CH4.tif fusion_method=[Linear Blending] regression_threshold=0.05 max/avg_displacement_threshold=1.00 absolute_displacement_threshold=1.00 compute_overlap subpixel_accuracy computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");

run("Stack to RGB");
saveAs("Tiff", outfile);
close();
