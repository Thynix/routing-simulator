set terminal png size 1200,800
set output "plot_hops_dist.png"

set xtics 5

set xlabel "Hops Taken"
set ylabel "Number of Occurences"

set style data linespoints
set style fill solid border -1

plot 'freenet-base-hops-policy-0' title 'Base Policy 0',\
     'freenet-base-hops-policy-1' title 'Base Policy 1',\
     'freenet-flat-hops-policy-0' title 'Flat Policy 0',\
     'freenet-flat-hops-policy-1' title 'Flat Policy 1'
