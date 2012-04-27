set terminal png size 1200,800
set output "walk_percentile.png"

set key left top

set title "Endpoints: |Network| Random Walks"

set xlabel "Hops"
set xrange [0:]

set ylabel "Percent Nodes Seen"

plot [0:50] 63.5825 title "Reference",\
    '90.dat' title "90th Percentile" with points,\
    '97.dat' title "97th Percentile" with points,\
    '99.dat' title "99th Percentile" with points,\
    'mean.dat' title "Mean" with points

