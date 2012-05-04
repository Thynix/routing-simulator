set terminal png size 1200,800
set output "probeDistribution.png"

set key left top

set xrange [0:]

set title "Degree Distribution"
set xlabel "Node number"
set ylabel "Cumulative Appearances"

plot 'occurenceDistribution/reference.dat' title "Reference" s cumul,\
     'occurenceDistribution/probe-50.dat' title "50 HTL" s cumul,\
     'occurenceDistribution/probe-45.dat' title "45 HTL" s cumul,\
     'occurenceDistribution/probe-40.dat' title "40 HTL" s cumul,\
     'occurenceDistribution/probe-35.dat' title "35 HTL" s cumul,\
     'occurenceDistribution/probe-30.dat' title "30 HTL" s cumul,\
     'occurenceDistribution/probe-25.dat' title "25 HTL" s cumul,\
     'occurenceDistribution/probe-10.dat' title "10 HTL" s cumul,\
     'occurenceDistribution/probe-5.dat' title "5 HTL" s cumul,\
     'occurenceDistribution/probe-2.dat' title "2 HTL" s cumul,\
     'occurenceDistribution/probe-1.dat' title "1 HTL" s cumul,\
     'occurenceDistribution/probe-0.dat' title "0 HTL" s cumul

