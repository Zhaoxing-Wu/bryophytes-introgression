#cd("/Users/zhaoxingwu/Desktop/claudia lab/2023 phylogenetics data analysis")


using PhyloNetworks, PhyloPlots, DataFrames, CSV, Statistics, Distributions, Random, DelimitedFiles, Combinatorics
using JLD2

#https://github.com/chaoszhang/Weighted-ASTRAL_data
#https://www.sciencedirect.com/science/article/pii/S0960982218311254

#choose(13, 4) = 715
#choose(12, 4)
function main()
    #tree_remove_species()
    run_snaq()
    
end

function run_snaq()
    cf = load_object("cf.jld2")
    @time net0 = snaq!(readTopology("tree_start_ord.txt"), readTableCF(cf), hmax=1, filename="net0_bucky", seed=123, runs=1)
end

# pick a tree as the starting tree for snaq
# convert the species name to order
function tree_remove_species()
    genetrees = readMultiTopology("./tree.txt"); # list of all trees

    tm = DataFrame(CSV.File("./group.csv"))
    taxonmap = Dict(row[:Species] => row[:Order] for row in eachrow(tm))

    file = open("tree_start_ord.txt", "a") #output

    for net in genetrees
        leaf = []
        for i in net.leaf
            push!(leaf, i.name)
        end

        flag = Dict(Pair.(tm[:, "Order"], false)) #whether the order exists in the tree
        net_merge = net #output tree
        leaf_del = [] #list of duplicated order

        cnt = 0 #count the number of leaf added
        for i in 1:length(leaf)
            if leaf[i] in tm[:,"Species"]
                if !flag[taxonmap[leaf[i]]] #found an order for output
                    flag[taxonmap[leaf[i]]] = true
                    net_merge.leaf[i].name = string(taxonmap[leaf[i]])
                    cnt += 1
                else #remove the leaf
                    push!(leaf_del, leaf[i])
                end
            end
        end

        for i in leaf_del #remove the leaf
            deleteleaf!(net_merge, i, keeporiginalroot=true)
        end
        
        if cnt == length(unique(tm[:, "Order"])) #if output tree contain all orders
            print(net_merge)
            write(file, writeTopology(net_merge))
            break
        end
    end
    close(file)
end

function read_cf()
    tm = DataFrame(CSV.File("./group.csv"))
    taxonmap = Dict(row[:Species] => row[:Order] for row in eachrow(tm))
    genetrees = readMultiTopology("./tree.txt");
    q,t = countquartetsintrees(genetrees, taxonmap, showprogressbar=true);
    save_object("q.jld2", q)
    save_object("t.jld2", t)
    df_wide = writeTableCF(q,t)
    df = df_wide[:,[:t1, :t2, :t3, :t4, :CF12_34, :CF13_24, :CF14_23]]
    save_object("cf.jld2", df)
end


# convert .nex file to a format that can be accepted by readMultiTopology()
function convert_nex_tree()
    f = open("./orig_data/SnAq_Bryophyte_Introgression/Astral_allbryos_G_Jul2021_GeneTrees_rooted.nex", "r")
    file = open("./tree.txt", "a")

    for line in readlines(f)
        if startswith(line, "\ttree")
            line_clean = lstrip(split(line, "]")[end])
            write(file, line_clean*"\n")
        end
    end
    close(f)
    close(file)
end

main()