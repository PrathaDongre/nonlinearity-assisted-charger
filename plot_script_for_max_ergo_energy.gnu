# Set output file type and name
set terminal pngcairo enhanced font 'Times new roman,16'
set output 'ergotropy_and_energy_vs_time_npnlin.png'

# Set title and labels
set title "Harmonic, non-linear(g_2=0.05) optimal charging time for different F."
set xlabel "F/g_1"
set ylabel "g_1 τ"

# Set plot style to points and crosses
set style data points
set style data linespoints

# Plot the data from max_ergotropy_time file
set key at graph 1.0, 0.9

# Plot the data from max_ergotropy_time file
plot 'max_ergotropy_time_g1_0.1_g2_0.05_F_diff.dat' using ($1/0.1):($3*0.1) with points pointtype 7 pointsize 2 title "Ergotropy",\
 'max_energy_time_g1_0.1_g2_0.05_F_diff.dat' using ($1/0.1):($3*0.1) with points pointtype 2 pointsize 2 title "Energy"