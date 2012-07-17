set terminal png size 1200,800
set output "plot_peer_dist.png"

set xrange [1:50]
set xtics 5

set xlabel "Reported Peers"
set ylabel "Nodes"

set style data linespoints
set style fill solid border -1

plot 'base-degree' title 'Initial',\
     'freenet-base-degree' title 'Freenet on Ideal',\
     'freenet-flat-degree' title 'Freenet on Flat'
