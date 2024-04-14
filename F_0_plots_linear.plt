# Set up the plot
set terminal pngcairo enhanced font 'Times new roman,20'
set output 'F_0_plots_linear.png'

# Set axis titles
set xlabel "Time"
set ylabel "E or {/Symbol e}"

set key at graph 1.0, 0.85


plot "energy_g1_0.1_g2_0.0_F_0.0.dat" u 1:2 w l lc rgb "#880D1E" lw 3 title "Energy,F=0.0",\
"ergotropy_g1_0.1_g2_0.0_F_0.0.dat" u 1:2 w l lc rgb "#8ea604" lw 3 title "Ergotropy,F=0.0"
