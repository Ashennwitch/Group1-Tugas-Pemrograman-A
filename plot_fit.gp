set terminal pngcairo size 800,600
set output 'fit_plot.png'
set title 'Data vs Linear Fit'
set xlabel 'X'
set ylabel 'Y'

# Tambahkan baris ini supaya Gnuplot tahu pakai koma sebagai pemisah
set datafile separator ","

# Plot data asli (data.csv) sebagai titik, dan hasil fitting sebagai garis
plot 'sample.csv'       using 1:2 with points pt 7 title 'Original', \
     'fit_output.csv' using 1:2 with lines  lw 2 title 'Linear Fit'
