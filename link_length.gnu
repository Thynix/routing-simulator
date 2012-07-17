set terminal png size 1200,800

set title 'Link Length Distribution'

set xlabel 'Link Length (delta location)'
set ylabel 'Percent nodes with this length or less'

set output "plot_link_length.png"
set logscale x

#As location is circular and [0,1), largest difference is 0.5.
plot [0.00001:0.5] [0:1] 'base-link' s cumul title 'Ideal',\
                         'flat-link' s cumul title 'Flat',\
                         'freenet-base-link' s cumul title 'Freenet on Ideal',\
                         'freenet-flat-link' s cumul title 'Freenet on Flat',\
                         '1407_links' s cumul title 'Measured'

