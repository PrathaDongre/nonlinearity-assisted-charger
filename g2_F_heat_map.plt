# Set up the plot
set terminal pngcairo enhanced font 'Times new roman,15'
set output 'test_max_g2_F_heatmap.png'

# Set data range
set xrange [0:7]  # Assuming 8 columns for g2 values
set yrange [0:7]  # Assuming 8 rows for F values

# Set axis labels
set xlabel "g2"
set ylabel "F"

# Set custom ticks for x-axis
set xtics ("1.5" 0, "1.3" 1, "1.1" 2, "0.9" 3, "0.7" 4, "0.5" 5, "0.3" 6, "0.1" 7)
set ytics ("0.5" 0, "0.4" 1, "0.3" 2, "0.25" 3, "0.2" 4, "0.15" 5, "0.1" 6, "0.05" 7)

# Set color palette
set palette defined (0 "#ffd930", 0.5 "#DD2D4A", 1 "#0d2f6e")

# Plot the max energy matrix
# plot 'max_ergotropy_matrix.dat' matrix with image title "Max ergotropy Matrix", 
plot 'steady_ergotropy_matrix.dat' matrix with image title "Steady Ergotropy Matrix"
 