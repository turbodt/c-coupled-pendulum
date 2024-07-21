set terminal pngcairo size 800,600 enhanced font 'Arial, 14'
set output './data/graph.png'

filedata = './data/out.txt'
n = system(sprintf('cat %s | wc -l', filedata))

set title "Line Graph of Two Columns"
set xlabel "X-axis Label"
set ylabel "Y-axis Label"

set grid

set style line 1 lt 1 lw 2 pt 7 ps 1.5

plot filedata using 1:2 with lines linestyle 1 title 'Data Line'

set output
