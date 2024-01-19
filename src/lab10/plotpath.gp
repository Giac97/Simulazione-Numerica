# Set the output file format and filename
set terminal pngcairo enhanced truecolor font "arial,12" size 1920,1080
set output "tsp.png"

# Set plot title and axis labels
set title "USA TSP"
set xlabel "Longitude"
set ylabel "Latitude"

# Define the data file and column mapping
datafile = "path.out"
x_column = 1  # Column number for X values
y_column = 2  # Column number for Y values

# Customize line style and data points
set style line 1 lc rgb "blue" lw 2  # Line color and width

# Plot the data
plot datafile using x_column:y_column with linespoints ls 1 title "Data Series"

