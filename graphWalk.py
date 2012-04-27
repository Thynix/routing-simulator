import re

#Takes as input the random walk percentile sample size. Other lines ignored.
#<hop> hops:
#Mean:[whitespace][numers].[numbers]
#<target>th percentile:[whitespace][numbers].[numbers]
#Writes to each percentile.dat with the <last seen hop> <value given>\n

input = open("partial-results.txt", "r")

hopRe = re.compile(r"^(\d+) hops:")
meanRe = re.compile(r"^Mean:\s+(\d+.\d+)")

mean = open("mean.dat", "w")

targetPercentiles = [ 90, 97, 99 ]
percentiles = []
for percentile in targetPercentiles:
    percentiles.append((re.compile(r"^{}th percentile:\s+(\d+.\d+)".format(percentile)),
                        open("{}.dat".format(percentile), "w")))

hop = 0
match = ""

def writeMatch(file, match):
    file.write("{} {}\n".format(hop, match[0]))

for line in input:
    match = hopRe.findall(line)
    if len(match) is 1:
        hop = match[0]
        continue
    
    match = meanRe.findall(line)
    if len(match) is 1:
        writeMatch(mean, match)
    
    for percentile in percentiles:
        match = percentile[0].findall(line)
        if len(match) is 1:
            writeMatch(percentile[1], match)

