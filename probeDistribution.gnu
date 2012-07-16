set terminal png size 1200,800
set output "probeDistribution.png"

set key left top

set xrange [0:]

#set title "Uniform Probe Distribution on Ideal Network"
#set title "Uniform Probe Distribution on Degree-Conforming Network"
#set title "Metropolis-Hastings Corrected Probe Distribution on Ideal Network
set title "Metropolis-Hastings Corrected Probe Distribution on Degree-Conforming Network"
#set title "Probe Distribution"

unset xtics
set ylabel "Cumulative Appearances"
unset ytics

plot 'occurenceDistribution/reference.dat' title "Reference" s cumul,\
     'occurenceDistribution/probe-10.dat' title "10 HTL" s cumul,\
     'occurenceDistribution/probe-9.dat' title "9 HTL" s cumul,\
     'occurenceDistribution/probe-8.dat' title "8 HTL" s cumul,\
     'occurenceDistribution/probe-7.dat' title "7 HTL" s cumul,\
     'occurenceDistribution/probe-6.dat' title "6 HTL" s cumul,\
     'occurenceDistribution/probe-5.dat' title "5 HTL" s cumul,\
     'occurenceDistribution/probe-4.dat' title "4 HTL" s cumul,\
     'occurenceDistribution/probe-3.dat' title "3 HTL" s cumul,\
     'occurenceDistribution/probe-2.dat' title "2 HTL" s cumul,\
     'occurenceDistribution/probe-1.dat' title "1 HTL" s cumul

