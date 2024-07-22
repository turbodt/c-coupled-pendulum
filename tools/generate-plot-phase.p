set terminal pngcairo size 800,600 enhanced font 'Arial, 14'
set output './data/phase.png'

filedata = './data/out.txt'
n = system(sprintf('cat %s | wc -l', filedata))

set key autotitle columnhead
set title "Line Graph of Two Columns"

set grid

set style line 1 lt 1 lw 2 pt 7 ps 1.5 lc rgb 'blue'
set style line 2 lt 1 lw 2 pt 7 ps 1.5 lc rgb 'red'

plot \
    filedata using 2:3 with lines linestyle 1 title 'theta_1', \
    filedata using 2:4 with lines linestyle 2 title 'theta_2'

set output
