set terminal png size 1200,800

set title 'Link Length Distribution (8,000 node network)'

set xlabel 'Link Length (delta location)'
set ylabel 'Percent nodes with this length or less'

set output "plot_link_length.png"

#As location is circular and [0,1), largest difference is 0.5.
plot [0:0.5] [0:1] 'initial' s cumul title 'Start: Kleinberg',\
                   'requests-80000' s cumul title '80k requests',\
                   'requests-160000' s cumul title '160k requests',\
                   'requests-320000' s cumul title '320k requests',\
                   'requests-640000' s cumul title '640k requests',\
                   '1407_links' s cumul title 'Measured'


set output "plot_link_length_log.png"
set logscale x

plot [0.000001:0.5] [0:1] 'initial' s cumul title 'Start: Kleinberg',\
                   'requests-80000' s cumul title '80k requests',\
                   'requests-160000' s cumul title '160k requests',\
                   'requests-320000' s cumul title '320k requests',\
                   'requests-640000' s cumul title '640k requests',\
                   '1407_links' s cumul title 'Measured'
