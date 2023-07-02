using PhyloNetworks, PhyloPlots, DataFrames, CSV, Statistics, Distributions, Random, DelimitedFiles, Combinatorics
using JLD2

# 47 taxa
t = ["Andreaeales", "Anthocerotales", "Archidiales", "Aulacomniales", "Bartramiales", "Blasiales", "Bruchiales", "Bryales", "Buxbaumiales", 
"Catoscopiales", "Dicranales", "Ditrichales", "Erpodiales", "Fossombroniales", "Funariales", "Gigaspermales", "Grimmiales", "Hedwigiales", 
"Hookeriales", "Hypnales", "Hypnodendrales", "Hypopterygiales", "Jungermanniales", "Marchantiales", "Metzgeriales", "Myliales", 
"Neohodgsoniales", "Oedipodiales", "Orthodontiales", "Orthotrichales", "Pallaviciniales", "Pelliales", "Perssoniales", "Pleurophascales", 
"Pleuroziales", "Polytrichales", "Porellales", "Pottiales", "Ptilidiales", "Ptychomniales", "Rhabdoweisiales", "Rhizogoniales", 
"Sphaerocarpales", "Sphagnales", "Splachnales", "Takakiales", "Tetraphidales"]

file_phy = open("data/phy.txt", "r")
df = DataFrame(CSV.File("data/group.csv"))
file_seq = open("hyde_data.txt", "a")
file_map = open("hyde_map.txt", "a")
for line in readlines(file_phy)
    species = split(line, "\t")[1]
    for row in eachrow(df)
        if row.Species == species 
            if row.Order in t
                write(file_map, species*"\t"*row.Order*"\n")
                write(file_seq, line*"\n")
            end
            break
        end
    end
end
close(file_seq)
close(file_map)
    
