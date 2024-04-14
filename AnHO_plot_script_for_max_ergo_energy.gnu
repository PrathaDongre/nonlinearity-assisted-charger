# Set output file type and name
set terminal pngcairo enhanced font 'Times new roman,16'
set output 'AnHO_ergotropy_and_energy_vs_time_nonlinear_g2_0.1.png'

# Set title and labels
set title "Anharmonic optimal charging time for different F."
set xlabel "F/g_1"
set ylabel "g_1 τ"

# Set plot style to points and crosses
set style data points
set style data linespoints

# Plot the data from max_ergotropy_time file
set key at graph 1.0, 0.9

plot 'max_ergotropy_time_g1_0.1_anHO_F_diff (copy).dat' using ($1/0.1):($3*0.1) with points pointtype 7 pointsize 2 title "Ergotropy",\
'max_energy_time_g1_0.1_anHO_F_diff.dat' using ($1/0.1):($3*0.1) with points pointtype 2 pointsize 2 title "Energy"