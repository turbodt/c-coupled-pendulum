set terminal pngcairo size 800,600 enhanced font 'Arial, 14'
set output './data/energy.png'

filedata = './data/out.txt'
n = system(sprintf('cat %s | wc -l', filedata))

set key autotitle columnhead
set title "System energy and potentials"

set grid

set style line 1 lt 1 lw 2 pt 7 ps 1.5 lc rgb 'blue'
set style line 2 lt 1 lw 2 pt 7 ps 1.5 lc rgb 'red'
set style line 3 lt 1 lw 2 pt 7 ps 1.5 lc rgb 'green'
set style line 4 lt 1 lw 2 pt 7 ps 1.5 lc rgb 'black'

plot \
    filedata using 2:19 with lines linestyle 1 title 'T', \
    filedata using 2:20 with lines linestyle 2 title 'U_k', \
    filedata using 2:21 with lines linestyle 3 title 'U_h', \
    filedata using 2:22 with lines linestyle 4 title 'E' \
;

set output
