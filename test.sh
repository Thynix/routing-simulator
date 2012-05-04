#/bin/bash
#Exit if any command exits with an error
set -o errexit

#Start from directory in which the script is located.
cd `dirname $0`
startingDir=`pwd`

###Configure###
#Location of the routing simulator
jar="$startingDir/target/routing_simulator-0.0.1-dev-jar-with-dependencies.jar"
#Path to pyProbe installation (for plotting gnuplot scripts)
gnuplotDir="$startingDir/../../stats"
#General network size
size=3000
#Ideal graph settings
ideal="--ideal --local 0 --remote 6"
#Degree graph settings
degree="--degree $startingDir/../../stats/peerDist.dat_1407 --force-size"
#Simulation settings
simulation="--probe 50"

###End configure###

#TODO: Some way to not repeat directory names, or generate them from simulation parameters. Perhaps Python would be better-suited?

java -jar "$jar" $degree --size "$size" --metropolis-hastings -D "degree_mh_$size/peerDist.dat" -L "degree_mh_$size/links_output" -O "degree_mh_$size/occurenceDistribution/" $simulation
java -jar "$jar" $degree --size "$size" -D "degree_uniform_$size/peerDist.dat" -L "degree_uniform_$size/links_output" -O "degree_uniform_$size/occurenceDistribution/" $simulation

java -jar "$jar" $ideal --size "$size" --metropolis-hastings -D "ideal_mh_$size/peerDist.dat" -L "ideal_mh_$size/links_output" -O "ideal_mh_$size/occurenceDistribution/" $simulation
java -jar "$jar" $ideal --size "$size" -D "ideal_uniform_$size/peerDist.dat" -L "ideal_uniform_$size/links_output" -O "ideal_uniform_$size/occurenceDistribution/" $simulation

for dir in "degree_mh_$size" "degree_uniform_$size" "ideal_mh_$size" "ideal_uniform_$size"
do
    cd "$dir"
    #Link length distribution
    gnuplot "$gnuplotDir"/link_length.gnu
    #Degree distribution
    gnuplot "$gnuplotDir"/peer_dist.gnu
    gnuplot "$startingDir"/probeDistribution.gnu
    cd "$startingDir"
done
