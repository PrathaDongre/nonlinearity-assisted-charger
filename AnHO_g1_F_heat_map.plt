# Set up the plot
set terminal pngcairo enhanced font 'Times new roman,15'
set output 'AnHO_max_g1_F_heatmap.png'

# Set data range
set xrange [0:7]  # Assuming 8 columns for g2 values
set yrange [0:7]  # Assuming 8 rows for F values

# Set axis labels
set xlabel "g_1"
set ylabel "F"


# Set custom ticks for x and y axes
set xtics ("0.93" 0, "0.81" 1, "0.69" 2, "0.57" 3, "0.45" 4, "0.33" 5, "0.21" 6, "0.09" 7)
set ytics ("1.1" 0, "0.9" 1, "0.7" 2, "0.5" 3, "0.3" 4, "0.1" 5, "0.09" 6, "0.07" 7)
# Set color palette
set palette defined (0 "#ffd930", 0.5 "#DD2D4A", 1 "#0d2f6e")

# Plot the max energy matrix
plot 'AnHO_max_ergotropy_matrix.dat' matrix with image title "Anharmonic max ergotropy Matrix", 
