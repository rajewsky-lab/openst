arglist = getArgument();
args = split(arglist,';');
var gridx = args[0];
var gridy = args[1];
var indir = args[2];
var outfile = args[3];
var filename_template = args[4];

run("Grid/Collection stitching", "type=[Grid: snake by rows] order=[Right & Down] grid_size_x=" + gridx + " grid_size_y=" + gridy + " tile_overlap=20 first_file_index_i=0 directory=[" + indir + "] file_names=" + filename_template + " fusion_method=[Linear Blending] regression_threshold=0.05 max/avg_displacement_threshold=1.00 absolute_displacement_threshold=1.00 compute_overlap subpixel_accuracy computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
run("RGB Color");
saveAs("Tiff", outfile);
close();