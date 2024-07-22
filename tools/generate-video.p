set terminal webp animate delay 10 size 512, 512;
unset border;
unset xtics;
unset ytics;
unset key;

set output './data/out.webp';

filedata = './data/out.txt'

n = system(sprintf('cat %s | wc -l', filedata))

set xrange [ -0.5:1.5]
set yrange [ -1.5:0.5]

max(a,b) = (a < b) ? b : a;

do for [j=1:n] {

  set title sprintf('t = %d', j)

  plot \
    filedata using 7:8 every ::j::j with linespoints lc "red", \
    filedata using 9:10 every ::j::j with linespoints lc "red", \
    filedata using 11:12 every ::j::j with linespoints lc "red", \
    filedata using 13:14 every ::j::j with linespoints lc "red", \
    filedata using 15:16 every ::j::j with linespoints lc "red", \
    filedata using 17:18 every ::j::j with linespoints lc "red" \
  ;
}

unset output
