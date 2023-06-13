using PhyloNetworks, PhyloPlots, DataFrames, CSV, Statistics, Distributions, Random, DelimitedFiles, Combinatorics
using JLD2

f = open("data/hyde_data.txt", "r")
df = DataFrame(CSV.File("data/group.csv"))
file = open("temp.txt", "a")
for line in readlines(f)
    species = split(line, "\t")[1]
    for row in eachrow(df)
        if row.Species == species
            write(file, species*"\t"*row.Order*"\n")
            break
        end
    end
end
close(f)
close(file)
    
