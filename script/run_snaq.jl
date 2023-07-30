#cd("/Users/zhaoxingwu/Desktop/claudia lab/2023 phylogenetics data analysis/bryophytes-introgression")


using PhyloNetworks, PhyloPlots, DataFrames, CSV, Statistics, Distributions, Random, DelimitedFiles, Combinatorics
using JLD2

#https://github.com/chaoszhang/Weighted-ASTRAL_data
#https://www.sciencedirect.com/science/article/pii/S0960982218311254

#choose(13, 4) = 715
#choose(12, 4)
function main()

    """#keep 47 taxa
    t = [#taxa in quartets with cf value 1 and 0
        "Andreaeobryales", "Calobryales", "Diphysciales", "Disceliales", "Encalyptales", "Flexitrichales", "Leiosporocerotales", "Timmiales", "Treubiales",
        #taxa in hornworts
        "Dendrocerotales", "Notothyladales", "Phymatocerotales", "Leiosporocerotales"]
    for i in t
        cf = filter(row -> !(row.t1 == i || row.t2 == i || row.t3 == i || row.t4 == i),  cf)
    end
    """

    """ #unique taxa
    temp = []
    for i in 1:size(cf)[1]
        push!(temp, cf[i, "t1"])
        push!(temp, cf[i, "t2"])
        push!(temp, cf[i, "t3"])
        push!(temp, cf[i, "t4"])
    end
    print(length(unique!(temp)))
    """
    #tree_remove_species()
    #run_snaq()
    snaq = "((((((((((Oedipodiales,Tetraphidales):0.0474189866674788,Polytrichales):0.11666835419656903,(((Takakiales,((((Pelliales,(Pallaviciniales,Fossombroniales):0.15176482021691703):0.051324947478111956,((Pleuroziales,Metzgeriales):0.8809347307610181,(((Myliales,Jungermanniales):1.4312258020893758,Perssoniales):0.5120616991814858,(Porellales,Ptilidiales):0.007509442728669249):0.6554168215068044):0.1637814133505386):1.0235510098539058,(Blasiales,((Marchantiales,Sphaerocarpales):0.044010709826757755,Neohodgsoniales):1.4431139041628718):0.8692234502142334):2.0106675329161967,Anthocerotales):1.4254545472741615):0.1645769401922199,Sphagnales):0.5843133916516633,Andreaeales):0.6418901625370256):0.4683347895138096,Buxbaumiales):0.6818145764795444,Gigaspermales):0.3595030560064854,Funariales):0.4728864283980148,(((Grimmiales,(((Rhabdoweisiales,((Pottiales,Ditrichales):0.14759290572486392,((Erpodiales,Bruchiales):0.5235225142216609)#H48:0.0074656727566994536::0.9248139135848961):0.5588610757593384):0.048648675205419975,(Dicranales,#H48:0.032825715845608276::0.07518608641510395):3.331340458707548):0.2135087048736687,Archidiales):0.059691338040607146):0.3959014350055196,Pleurophascales):0.2847956782061885,Catoscopiales):0.48460063445423585):1.2245905859534711,(Bryales,((Rhizogoniales,(((Aulacomniales,Orthodontiales):0.04556777293890821,Hypnodendrales):0.0011262255415632537,((Hypopterygiales,(Hookeriales,Hypnales):0.045278584169739136):0.7429364702696049,Ptychomniales):0.5624696028657126):0.14064609375443174):0.00024867380288660343,Orthotrichales):0.5543727265862737):0.40720992080206314):6.65414557359528e-5,Splachnales):0.07783438026976097,Bartramiales,Hedwigiales);"
    plot(readTopology(snaq), :R)#, style=:fulltree, showGamma=false, useEdgeLength=true, showEdgeLength=true)
end

function run_snaq()
    cf = load_object("./data/cf_removed.jld2")
    @time net0 = snaq!(readTopology("./data/tree_start_ord.txt"), readTableCF(cf), hmax=1, filename="net0_bucky", seed=123, runs=1,
    ftolRel = 1e-4, ftolAbs = 1e-4, xtolRel = 1e-2, xtolAbs = 1e-2, liktolAbs = 1e-4)
end


# convert the species name to order
function tree_remove_species()
    genetrees = readMultiTopology("./data/tree.txt"); # list of all trees

    tm = DataFrame(CSV.File("./data/group.csv"))
    taxonmap = Dict(row[:Species] => row[:Order] for row in eachrow(tm))

    #tree_start_ord
    file = open("tree_ord.txt", "a") #output
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
                if taxonmap[leaf[i]] == "Dendrocerotales" || taxonmap[leaf[i]] == "Notothyladales" || taxonmap[leaf[i]] == "Phymatocerotales" || taxonmap[leaf[i]] == "Leiosporocerotales"    
                    #remove other leaf in hornwort order
                    push!(leaf_del, leaf[i])
                elseif !flag[taxonmap[leaf[i]]] #found an order for output
                    flag[taxonmap[leaf[i]]] = true
                    net_merge.leaf[i].name = string(taxonmap[leaf[i]])
                    cnt += 1
                else #remove the leaf of duplicated order
                    push!(leaf_del, leaf[i])
                end  
            end
        end

        for i in leaf_del #remove the leaf
            deleteleaf!(net_merge, i, keeporiginalroot=true)
        end
        
        """
        # pick a tree as the starting tree for snaq
        if cnt == length(unique(tm[:, "Order"])) #if output tree contain all orders
            print(net_merge)
            write(file, writeTopology(net_merge))
            break
        end
        """
        write(file, writeTopology(net_merge)*"\n")
    end
    close(file)
end

function read_cf()
    tm = DataFrame(CSV.File("./data/group.csv"))
    taxonmap = Dict(row[:Species] => row[:Order] for row in eachrow(tm))
    genetrees = readMultiTopology("./data/tree.txt");
    q,t = countquartetsintrees(genetrees, taxonmap, showprogressbar=true);
    save_object("./data/q.jld2", q)
    save_object("./data/t.jld2", t)
    df_wide = writeTableCF(q,t)
    df = df_wide[:,[:t1, :t2, :t3, :t4, :CF12_34, :CF13_24, :CF14_23]]
    save_object("./data/cf.jld2", df)
end


# convert .nex file to a format that can be accepted by readMultiTopology()
function convert_nex_tree()
    f = open("./data/orig_data/SnAq_Bryophyte_Introgression/Astral_allbryos_G_Jul2021_GeneTrees_rooted.nex", "r")
    file = open("./data/tree.txt", "a")

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