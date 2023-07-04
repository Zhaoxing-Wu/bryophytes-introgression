 # author: Marianne
 # from: https://github.com/mbjorner/hybrid-detection-comparison/blob/main/scripts/SummaryTables/analyzeMethodOutputFile.jl
 # following are functions to ascertain the ghost lineages at root of the tree,
 # and remove them to see how well functions can discern ghost lineages vs. discern hybrids 
 # only within a triplet or quartet

 # this script takes a single filename and creates an output file that contains
 # true negatives, true positives, false positives, false negatives at alpha levels of 0.05 and the bonferroni correction

 # will also contain a column for ghost lineages

 using DataFrames, CSV, PhyloNetworks, Statistics, Distributions
 
 #filepath = file of hyde output to analyze
 #filename = desired filename for output to save
 #netname = filepath to proposed species tree or network
 #savefilepath = desired path to save output
 analyzeHyDe_DStat("./output/hyde-out.txt", 
 "temp", "./data/tree.txt", "./output", "Anthocerotales") 

 #include("/Users/bjorner/GitHub/hybrid-detection-comparison/scripts/HelperFunctions/calcD.jl")
 
 #net_directory = "/Users/bjorner/GitHub/ms_gene_trees/nets_newick/"

 #=
inspired from evobiR::CalcD but requires package "seqinr" only
used in [Blair & Ané 2020](https://doi.org/10.1093/sysbio/syz056)
=#

using BioSymbols
using PhyloNetworks
using Random
using Statistics
using Distributions
using DataFrames
using CSV

"""
    calcPatternCounts(filename; taxa=["P1","P2","P3","O"], ambig=:D)

Take a single 4-taxon alignment and calculate 3 numbers of sites:
total, ABBA sites, and BABA sites. By default, sites with any ambiguity
are dropped (`ambig=:D`) before calculating these 3 numbers.
Alternatively, they could be resolved (`ambig=:R`). Otherwise,
ambiguous sites are just left alone (e.g. ignored: `ambig=:I`).
In that case for instance, `N--N` would count as an ABBA pattern.

warning: bug in `PhyloNetworks.readFastaToArray` with some ambiguous sites,
due to a bug in BioSequences, see
[FASTA.Reader](https://github.com/BioJulia/BioSequences.jl/issues/64).
A quick & dirty trick is used here.
Fix this bug in `PhyloNetworks.readFastaToArray`: use
sequence(DNASequence, record) instead of sequence(record)
"""
function calcPatternCounts(file::String; taxa=[1,2,3,10], ambig=:D)
  # read alignment and make it a matrix
  species, alg = PhyloNetworks.readFastaToArray(file)
  ## bug: due to [FASTA.Reader](https://github.com/BioJulia/BioSequences.jl/issues/64)
  ## quick fix, as suggested in the github issue:
  for i in 1:length(alg)
    if !(eltype(alg[i]) <: DNA)
        alg[i] = DNASequence(convert(String, alg[i]))
    end
    print(length(alg[i]))
  end
  # alg = PhyloNetworks.readfastatodna(file) with next registered version
  # perhaps do better using BioSequences or BioAlignments directly:
  # r = FASTA.Reader(open("trial0_concat.fasta"))

  # check that there are 4 taxa, and all have same alignment length
  # sometimes, 1 or more taxa might be missing a gene of interest
  goodalign = length(alg) == 4
  ind = indexin(taxa, species) # species[ind] in correct order to look at ABBA and BABA
  nsites = length(alg[1])

  if goodalign
    goodalign = [length(alg[i]) for i in 2:4] == [nsites, nsites, nsites]
  end
  goodalign || return (missing,missing,missing)

  ## deal with ambiguity in sequence data
  # R: A or G
  # Y: C or T
  # S: G or C
  # W: A or T
  # K: G or T
  # M: A or C

  if ambig in [:D, :R]
    # D: drop all ambiguous sites
    # R: resolve, but first drop any site having a 3- or 4-ambiguity
    target = (BioSymbols.ACGT..., DNA_R, DNA_Y, DNA_S, DNA_W, DNA_K, DNA_M)
    for j in nsites:-1:1
      keepj = (ambig==:D ?
               all(s -> BioSymbols.iscertain(s[j]), alg) :
               all(s -> s[j] in target, alg) )
      if keepj continue; end
      for s in alg
        deleteat!(s,j)
      end
    end
  end
  nsites = length(alg[1]) # possibly less than before
  ## resolve remaining ambigous sites randomly
  if ambig == :R
    function resolver(x)
      if iscertain(x) return x; end
      if x==DNA_R return rand([DNA_A, DNA_G]); end
      if x==DNA_Y return rand([DNA_C, DNA_T]); end
      if x==DNA_S return rand([DNA_G, DNA_C]); end
      if x==DNA_W return rand([DNA_A, DNA_T]); end
      if x==DNA_K return rand([DNA_G, DNA_T]); end
      if x==DNA_M return rand([DNA_A, DNA_C]);
      else error("weird base: $x"); end
    end
    for a in alg
      for j in 1:nsites
        isambiguous(a[j]) || continue
        a[j] = resolver(a[j])
      end
    end
  end
  # now calculate the number of sites with pattern ABBA or BABA
  abba = 0
  baba = 0
  nsites>0 || return (0,0,0)
  pattern = Vector{DNA}(undef, 4) # to avoid garbage collection
  for j in 1:nsites
    # pattern = [alg[ind[i]][j] for i in 1:4]
    for i in 1:4 pattern[i] = alg[ind[i]][j]; end
    length(unique(pattern)) == 2 || continue   # focus on biallelic sites
    pattern[1] != pattern[2] || continue   # focus on xy.. patterns
    pattern[3] != pattern[4] || continue   # focus on ..xy patterns 
    if pattern[3] == pattern[1]
      baba += 1
    else # pattern[3] == pattern[2]
      abba += 1
    end
  end
  return nsites,abba,baba
end

"""
    calcD(x, y)

Take a vector `x` of ABBA counts and a vector `y` of BABA counts, assuming that
in each list: 1 count corresponds to 1 gene. Then calculate the D statistics for
the concatenated data, that is: [sum(ABBA) - sum(BABA)] / [sum(ABBA) + sum(BABA)]
"""
function calcD(x, y)
  sx = sum(x); sy = sum(y)
  return (sx - sy)/(sx + sy)
end

"""
     calcDztest_concatenated(x,y)
     calcDsd_concatenated(x, y, nbootstrap=500)

Take the outcome of `calcPatternCounts` from a concatenated alignment
(except for the total number of non-ignored sites):
number of ABBA sites `x` and number of BABA sites `y`, and test whether
the data is consistent with mean D=0.

Output: D statistic, its standard deviation, and test p-value (2-tailed).

Assumptions:

1. independent sites in the single (concatenated) alignment
  from which `x,y` come from.
2. the number of "informative" sites for the ABBA-BABA test, that is, the
  number of ABBA + BABA sites is considered fixed. In other words, the test
  conditions on the number of informative sites.

Note that D = (x-y)/(x+y) becomes phat - (1-phat) = 2*(phat -0.5)
where phat = estimated proportion of ABBA sites among informative sites,
that is, among ABBA & BABA sites.
Testing "mean D=0" is equivalent to testing p=0.5, which amounts to doing
is simple binomial test, or chi-square goodness-of-fit test if all sites
are assumed independent.

The first version does a z test (equivalent to chi-square goodness-of-fit test).

The second version calculates the standard deviation of the D statistics
using bootstrap, re-sampling the ABBA-BABA sites (*not* resampling all sites:
see assumption 2 above). Then calculates z = D/sd and
pvalue = two-tailed probability beyond z under N(0,1).
This second version could be modified to allow for linkage between sites
(in which case the z-test becomes invalid).
"""
function calcDztest_concatenated(x,y)
  n = x+y
  D = (x-y)/n
  phat = x/n
  z = (phat-0.5)*sqrt(n)*2 #  = (phat - p0) / sqrt(p0 * (1-p0) * 1/n) when p0=0.5
  pvalue = ccdf(Normal(), abs(z))*2
  return D, z, pvalue
end

function calcDsd_concatenated(x, y, nbootstrap=500)
  n = x+y
  D = (x-y)/n
  p_abba = x/n
  d = Distributions.Binomial(n, p_abba)
  bx = rand(d, nbootstrap)
  bD = (2*bx .- n) ./ n # y=n-x so x-y = x-(n-x) = 2x-n
  z = D / Statistics.std(bD)
  pvalue = ccdf(Normal(), abs(z))*2
  return D, z, pvalue
end

"""
    calcDsd_byloci(x, y, nbootstrap=500)

Take a list of ABBA counts and a list of BABA counts, like `calcD`.
Calculate the standard deviation of the D statistics using bootstrap,
re-sampling loci. This is assuming independent loci
(but not independent sites within loci).

Default is 500 bootstrap replicates.
"""
function calcDsd_byloci(x, y, nbootstrap=500)
  nloci = length(x)
  bD = Vector{Float64}(undef, nbootstrap)
  for b in 1:nbootstrap
    bsx=0; bsy=0
    bloci = rand(1:nloci, nloci) # sample loci with replacement
    for locus in bloci
      bsx += x[locus]
      bsy += y[locus]
    end
    bD[b] = (bsx - bsy)/(bsx + bsy)
  end
  return Statistics.std(bD)
end
 
 function loadNet(netname::AbstractString)
     #net = readTopology(string(net_directory, netname))
     net = readTopology(netname)
     return net
 end
 
 function calculate_WC(testReturnsHybrid, isTrueHybridID, isHybridInTriple)
     return testReturnsHybrid && !isTrueHybridID && isHybridInTriple
 end
 
 function calcDsd_concatenated(abba, baba, bbaa, outgroup_position, nbootstrap=500)
     if outgroup_position == 2
         #AABB = ABBA, ABBA = BABA
         baba = abba
         abba = bbaa
     elseif outgroup_position == 3
         # AABB = ABBA, ABAB = BABA
         abba = bbaa
         #baba = baba
     elseif outgroup_position == 1
         #baba = baba
         #abba = abba
     end
 
     n = abba+baba
     D = (abba-baba)/n
     p_abba = abba/n
     d = Distributions.Binomial(n, p_abba)
     bx = rand(d, nbootstrap)
     bD = (2*bx .- n) ./ n # y=n-x so x-y = x-(n-x) = 2x-n
     z = D / Statistics.std(bD)
     pvalue = ccdf(Normal(), abs(z))*2
     return D, z, pvalue
   end
 
 function calcDpsd_concatenated(abba, abab, bbaa, outgroup_position, nbootstrap=500)
     # if outgroup position is 3, proceed normally. else, change whether we use abba or aabb. positions 1 and 2 give symmetric results
     if outgroup_position == 2
         #AABB = ABBA, ABBA = BABA
         tempBABA = abab
         abab = abba
         abba = bbaa
         bbaa = tempBABA
     elseif outgroup_position == 3
         # AABB = ABBA, ABAB = BABA
         tempABBA = abba
         abba = bbaa
         bbaa = tempABBA
         #baba = baba
     elseif outgroup_position == 1
         #baba = baba
         #abba = abba
     end
     n = abba+abab+bbaa
     Dp = (abba-abab)/n
     p_abba = abba/n
     d = Distributions.Binomial(n, p_abba)
     bx = rand(d, nbootstrap)
     bD = (2*bx .- n) ./ n # y=n-x so x-y = x-(n-x) = 2x-n
     z = Dp / Statistics.std(bD)
     pvalue = ccdf(Normal(), abs(z))*2
     return Dp, z, pvalue
 end
 
 function calculate_FP_FN_TP_TN(testReturnsHybrid, isTrueHybridID)
     true_positive = testReturnsHybrid && isTrueHybridID
     false_negative = !testReturnsHybrid && isTrueHybridID
     true_negative = !testReturnsHybrid && !isTrueHybridID
     false_positive = testReturnsHybrid && !isTrueHybridID
 
     return false_positive, false_negative, true_positive, true_negative
 end
 
 #=
 returns whether tip subsets contain a hybrid relationship by removing nodes from a network
 =#
 function isHybrid!(net, taxa)
     for l in tipLabels(net)
         if l ∉ taxa
             PhyloNetworks.deleteleaf!(net,l,keeporiginalroot=true, simplify=true);
             #@show net
         end 
     end 
 
     if net.numHybrids > 0
         return true
     else
         return false
     end
 end
 
 #=
 Given an input net, returns a Set of tip names that arose from a hybridization event
 =#
 function getHybrids(net::HybridNetwork)
     hybrids = Set()
     for i in 1:net.numHybrids
         # hardwiredCluster is invariant to the hybrid edge chosen, returns identical results no matter the edge
         hybridclade = hardwiredCluster(net.hybrid[i].edge[1], tipLabels(net));
         trueHybrids = tipLabels(net)[hybridclade]
         for hyb in trueHybrids
             push!(hybrids, hyb)
         end
     end
     return hybrids
 end
 
 #=
 Returns if the specified hybrid is the actual (or one of the) hybrid(s) in the network
 =#
 function isTrueHybrid(net::HybridNetwork, HybridStringName::String)
     return issubset([HybridStringName], getHybrids(net))
 end
 
 #=
 Returns the depth of the hybrid_minor_edge_parent
 
 get height of hybrid node
 get height of tip
 tip - hybridnode
 =#
 function calculateLengthToHybridEdge(net, hybridLabel)
     directEdges!(net); preorder!(net);
     heights = PhyloNetworks.getHeights(net)
     tipHeight = findmax(heights)[1]
     heightVector = [node.number => (heights[i], node.name) for (i,node) in enumerate(net.nodes_changed)]
     hybridnumber = net.hybrid[hybridLabel].number
     hybridHeight = Dict(heightVector)[hybridnumber][1]
 
     return tipHeight - hybridHeight
 end
 
 #= HyDe output consists of the following headers describing site pattern frequencies.
 AAAA	AAAB	AABA	AABB	AABC	ABAA	ABAB	ABAC	ABBA	BAAA	ABBC	CABC	BACA	BCAA	ABCD
 To replace analysis of site patterns using classical D3, site patterns can be converted to pairwise distances
 where d1-3, d2-3 can be computed by:
 
 d1-3
 sum over -> AABA AABB AABC CABC BACA BCAA ABCD ABBA BAAA ABBC 
     AAAB		AABB	AABC	ABAA		ABAC	ABBA		ABBC	CABC		BCAA	ABCD
 d2-3
 sum over -> AABA AABB AABC CABC BACA BCAA ABCD ABAB ABAA ABAC 
     AAAB	AABA		AABC		ABAB	ABAC	ABBA		ABBC	CABC	BACA		ABCD
 
 
 where AABA AABB AABC CABC BACA BCAA ABCD are equal between d1-3 and d2-3
 =#
 function calcD3fromHyDe(AAAB, AABA,	AABB,	AABC,	ABAA,	ABAB,	ABAC,	ABBA,	BAAA,	ABBC,	CABC,	BACA,	BCAA,	ABCD, relative_outgroup, nbootstrap=500)
     # if relative outgroup is in position 3, proceed normally
     d13, d23, d13, d12 = 0, 0, 0, 0
 
     if relative_outgroup == 3
         d13 = AAAB + AABB + AABC + ABAA + ABAC + ABBA + ABBC + CABC + BCAA + ABCD
         d23 = AAAB + AABA + AABC + ABAB + ABAC + ABBA + ABBC + CABC + BACA + ABCD
     elseif relative_outgroup == 2 #: d13 = d12, d23 = d23
         d12 = AABA + AABB + AABC + ABAA + ABAB + ABAC +	CABC + BACA + BCAA + ABCD
         d23 = AAAB + AABA + AABC + ABAB + ABAC + ABBA + ABBC + CABC + BACA + ABCD
     elseif relative_outgroup == 1
         d13 = AAAB + AABB + AABC + ABAA + ABAC + ABBA + ABBC + CABC + BCAA + ABCD
         d12 = AABA + AABB + AABC + ABAA + ABAB + ABAC +	CABC + BACA + BCAA + ABCD
     end
 
     # if relative outgroup is in position 1, compute D3 with d12 (== d23) and d13.
     if relative_outgroup == 1      
         d23 = d12 # so when +, there is an excess between 1 & 3, -, there is excess between 1 & 2
     elseif relative_outgroup == 2 # compute D3 with d12 and d23
         d13 = d12 # so when +, there is an excess between 1 & 2, -, there is excess between 2 & 3
     end
 
     n = d13 + d23
     D3 = (d13 - d23) / n
     p_d13 = d13 / n
     d = Distributions.Binomial(n, p_d13)
     bx = rand(d, nbootstrap)
     bD3 = (2*bx .- n) ./ n # y=n-x so x-y = x-(n-x) = 2x-n
     z = D3 / Statistics.std(bD3)
     pvalue = ccdf(Normal(), abs(z))*2
 
     return D3, z, pvalue
 end
 
 
 function makeUltrametric!(net::HybridNetwork)
     ## Base case: does not change the network as it is already ultrametric
     # if isUltrametric(net)
     #    return
     # end
 
     ## network is not ultrametric
     preorder!(net) # PhyloNetworks.directEdges!(net) ?? should we also direct the edges or is this unneeded
     heights = PhyloNetworks.getHeights(net)
     max_tip_height = findmax(heights)[1]
     heights_dict = Dict([node.name => (heights[i], node.number, node) for (i,node) in enumerate(net.nodes_changed)])
     labels = tipLabels(net)
     for tip_name in labels
         tip_height = heights_dict[tip_name][1]
         node = heights_dict[tip_name][3]
         if tip_height != max_tip_height # should there be floating point error accounted for?
             # change the length of the external edge that the tip is attached to to make it same length
             length_difference = max_tip_height - tip_height
             edgeToChange = PhyloNetworks.getMajorParentEdge(node)
             edgeToChange.length = edgeToChange.length + length_difference
         end
     end
 
     return net
 
 end
 
 
 #=
 Returns Tip:Root / Hybrid:Root ratio assuming ultrametric tree with all nodes
 =#
 function calculateTipToHybridRatio(net, hybridLabel)
     directEdges!(net); preorder!(net);
     heights = PhyloNetworks.getHeights(net)
     tipHeight = findmax(heights)[1]
     heightVector = [node.number => (heights[i], node.name) for (i,node) in enumerate(net.nodes_changed)]
     hybridnumber = net.hybrid[hybridLabel].number
     hybridHeight = Dict(heightVector)[hybridnumber][1]
 
     return tipHeight / hybridHeight
 end
 
 #=
 Returns the true mixing parameter of the minor edge of the given hybridlabel
 =#
 function getTrueMixingParameter(net, hybridLabel)
     minor_edge = PhyloNetworks.getMinorParentEdge(net.hybrid[hybridLabel])
     return minor_edge.gamma
 end
 
 # this indicates a ghost lineage at the root 
 function isGhostHybridAtRoot(net)
     rootEdges = net.node[net.root].edge
     for i in 1:net.numHybrids
         for n in net.hybrid[i].edge
             for r in rootEdges
                 if r.number == n.number
                     return true
                 end
             end
         end
     end
     return false
 end
 
 function removeHybridAtRoot(net)
     rootEdges = net.node[net.root].edge
     print("root is", rootEdges)
     for i in 1:net.numHybrids
         for n in net.hybrid[i].edge
             print("n.number is ", n.number)
             for r in rootEdges
                 if r.number == n.number
                     print("removing hybrid edge ", net.hybrid[i].edge)
                     PhyloNetworks.deletehybridedge!(net, n)
                 end
             end
         end
     end
 end
 
 function isGhostHybrid(net, hybridlabel)
     # checks to see if the hybrid clades of the network are downstream
     # of the hybrid edge
     PhyloNetworks.directEdges!(net)
     hybridclade = hardwiredCluster(net.hybrid[hybridlabel].edge[1], tipLabels(net))
     hybridtips = tipLabels(net)[hybridclade]
     print("\n hybrid tips are: ", hybridtips)
     # step 2 is find the clades downstream of the node / downstream of the hybrid
     # parent of minor hybrid edge:
     parentnode = PhyloNetworks.getParent(net.hybrid[hybridlabel].edge[1]) 
     print(parentnode)
     
     # children = PhyloNetworks.getChildren(parentnode)
     # majortree = PhyloNetworks.majorTree(net)
     parentparentnode.getChildren
 
     downstreamclades = PhyloNetworks.hardwiredCluster(parentnode.edge[3], tipLabels(net))
     downstreamtips = tipLabels(net)[downstreamclades]
     print("\n downstreamtips are", downstreamtips)
     if issubset(hybridtips, downstreamtips) 
         return true
     else
         # print("\n hybrid tips are ", hybridtips)
         # print("\n downstream tips are ", downstreamtips)
         return false
     end
 end
 
 #=
 The D3 statistic intends to be tested on events where the sister taxa are distinct from the 
 outgroup. D3 will behave differently (give false positives) if the underlying species tree
 is different than expected. Here, we test for a topology matching ((tax1,tax2),out) where
 both taxon 1 and taxon 2 (labels of leaves) are sister with the common outgroup "out" by
 ensuring the distance in the network's major tree from tax1->out == tax2->out
 =#
 function isSisterSisterOutgroup(tax1,tax2,out,net::HybridNetwork)
 
     majorNet = majorTree(net)
     majorNet = makeUltrametric!(majorNet)
     for l in tipLabels(net)
         if l ∉ [tax1, tax2, out]
             PhyloNetworks.deleteleaf!(majorNet,l,keeporiginalroot=true, simplify=true);
             #@show net
         end 
     end 
     taxonPairwiseMatrix = PhyloNetworks.pairwiseTaxonDistanceMatrix(majorNet)
     tipOrder = tipLabels(majorNet)
 
 
     tax1TipOrder = findall(x->x==tax1, tipOrder)[1]
     tax2TipOrder = findall(x->x==tax2, tipOrder)[1]
     outTipOrder = findall(x->x==out, tipOrder)[1]
 
     distTax1ToOut = taxonPairwiseMatrix[tax1TipOrder*outTipOrder]
     distTax2ToOut = taxonPairwiseMatrix[tax2TipOrder*outTipOrder]
    
     return abs(distTax1ToOut - distTax2ToOut) < 0.01
 end
 
 
 #=
 get length of minor hybrid edge of specified hybrid
 =#
 function getLengthMinorEdge(net, hybridlabel)
     minoredge = PhyloNetworks.getMinorParentEdge(net.hybrid[hybridlabel])
     return PhyloNetworks.getlengths([minoredge])
 end
 
 function isSignificant(column, alpha::Float64)
     return column .< alpha
 end
 
 function isSignificantBonferroni(column, alpha::Float64)
     bonferroni = alpha / size(column, 1)
     return column .< bonferroni
 end
 
 function analyzeHyDe_DStat(filepath, filename, netname, savefilepath) 
     HyDeOut = DataFrame(CSV.File(filepath))
     net = loadNet(netname)
     outFile = string(savefilepath, filename, "_analyzed.csv")
 
     sig05 = isSignificant(HyDeOut[!, :Pvalue], 0.05)
     sig_bonferroni = isSignificantBonferroni(HyDeOut[!, :Pvalue], 0.05)
     outgroup = split(filename,"h")[1]
     outgroup = split(outgroup,"n")[1]
     setsOfTriplets = size(HyDeOut,1)
 
     column_names = [:P1, :Hybrid, :P2, :HyDe_pvalue, :FP, :FN, :TP, :TN, :WC, :FP_bon, :FN_bon, :TP_bon, :TN_bon, :WC_bon, 
                     :hybridTripletExpected, :proposedHybridCorrect, :isGhostHybrid, :isGhostAtRoot, 
                     :isMultiHybrid, :lengthToHybridEdge, :tipToHybridRatio, 
                     :D, :D_z, :D_pvalue, :Dp, :Dp_z, :Dp_pvalue, :D3, :D3_z, :D3_pvalue, 
                     :FP_D, :FN_D, :TP_D, :TN_D, :FP_Dbon, :FN_Dbon, :TP_Dbon, :TN_Dbon, 
                     :FP_Dp, :FN_Dp, :TP_Dp, :TN_Dp, :FP_Dpbon, :FN_Dpbon, :TP_Dpbon, :TN_Dpbon, 
                     :FP_D3, :FN_D3, :TP_D3, :TN_D3, :FP_D3bon, :FN_D3bon, :TP_D3bon, :TN_D3bon, 
                     :true_gamma, :proposed_gamma, :isDtestTopology]
 
     # add all names to dataframe as empty columns
     df_results = DataFrame(column_names .=> Ref([]))
 
     for row in 1:setsOfTriplets
 
         net = loadNet(netname);
         P1 = string(HyDeOut[row, :P1]);
         Hybrid = string(HyDeOut[row, :Hybrid]);
         P2 = string(HyDeOut[row, :P2]);
         proposed_gamma = HyDeOut[row, :Gamma]
         true_gamma = -1
         triplet = [P1,Hybrid,P2];
         quad_includingoutgroup = [P1,Hybrid,P2,outgroup]
 
         hybridTripletExpected = isHybrid!(net, triplet)
         proposedHybridCorrect = false
         isGhostHybridBool = false
         isMultiHybrid = false
         isGhostHybridAtRootBool = false
 
         lengthToHybridEdge, tipToHybridRatio = -1, -1
 
         for i in 1:net.numHybrids
             if i == 2
                 isMultiHybrid = true
             end
             lengthToHybridEdge = calculateLengthToHybridEdge(net, i)
             tipToHybridRatio = calculateTipToHybridRatio(net, i)
             true_gamma = getTrueMixingParameter(net, i)
             proposedHybridCorrect = isTrueHybrid(net, Hybrid)
             isGhostHybridBool = isGhostHybrid(net, i) # is there a bug here?
             isGhostHybridAtRootBool = isGhostHybridAtRoot(net) # this means that there is a hybrid present not ""between"" tested taxa
         end
 
         FP, FN, TP, TN = calculate_FP_FN_TP_TN(sig05[row], proposedHybridCorrect)
         FP_bon, FN_bon, TP_bon, TN_bon = calculate_FP_FN_TP_TN(sig_bonferroni[row], proposedHybridCorrect)
         # WC = true if hyde identifies a hybrid in the triple as expected, but it's NOT the correct hybrid
 
         WC = calculate_WC(sig05[row], proposedHybridCorrect, hybridTripletExpected)
         WC_bon = calculate_WC(sig_bonferroni[row], proposedHybridCorrect, hybridTripletExpected)
 
         ABBA_site = HyDeOut[row,:ABBA]
         ABAB_site = HyDeOut[row,:ABAB]
         BBAA_site = HyDeOut[row,:AABB]
 
         AAAB = HyDeOut[row,:AAAB]
         AABA = HyDeOut[row,:AABA]
         AABB = HyDeOut[row,:AABB]	
         AABC = HyDeOut[row,:AABC]	
         ABAA = HyDeOut[row,:ABAA]	
         ABAB = HyDeOut[row,:ABAB]	
         ABAC = HyDeOut[row,:ABAC]	
         ABBA = HyDeOut[row,:ABBA]	
         BAAA = HyDeOut[row,:BAAA]	
         ABBC = HyDeOut[row,:ABBC]	
         CABC = HyDeOut[row,:CABC]	
         BACA = HyDeOut[row,:BACA]	
         BCAA = HyDeOut[row,:BCAA]	
         ABCD = HyDeOut[row,:ABCD]
 
         # D-statistics should only be calculated where the correct topology is followed. 
         isCorrectDTopology_pos3 = isSisterSisterOutgroup(P1, Hybrid, P2, net)
         isCorrectDTopology_pos1 = isSisterSisterOutgroup(P2, Hybrid, P1, net)
 
         outgroup_position = 0
         if isCorrectDTopology_pos1 
             outgroup_position = 1
         elseif isCorrectDTopology_pos3
             outgroup_position = 3
             # triplet is in the correct order
         else
             outgroup_position = 2
         end
 
         D, D_z, D_pvalue = calcDsd_concatenated(ABBA_site, ABAB_site, BBAA_site, outgroup_position)
         Dp, Dp_z, Dp_pvalue = calcDpsd_concatenated(ABBA_site, ABAB_site, BBAA_site, outgroup_position)
         D3, D3_z, D3_pvalue = calcD3fromHyDe(AAAB, AABA,	AABB,	AABC,	ABAA,	ABAB,	ABAC,	ABBA,	BAAA,	ABBC,	CABC,	BACA,	BCAA,	ABCD, outgroup_position)
         HyDe_pval = HyDeOut[row,:Pvalue]
         # bootstrap option to calculate D, z-score, and p-value
         #=
         FP_D, FN_D, TP_D, TN_D = calculateD_FP_FN_TP_TN(D, hybridTripletExpected, 0.05 > D_pvalue, net)
         FP_Dbon, FN_Dbon, TP_Dbon, TN_Dbon = calculateD_FP_FN_TP_TN(D, hybridTripletExpected, (0.05 * 3 / setsOfTriplets) > D_pvalue, net)
 
         FP_Dp, FN_Dp, TP_Dp, TN_Dp = calculateD_FP_FN_TP_TN(Dp, hybridTripletExpected, 0.05 > Dp_pvalue, net)
         FP_Dpbon, FN_Dpbon, TP_Dpbon, TN_Dpbon = calculateD_FP_FN_TP_TN(Dp, hybridTripletExpected, (0.05 * 3 / setsOfTriplets) > Dp_pvalue, net)
 
         FP_D3, FN_D3, TP_D3, TN_D3 = calculateD_FP_FN_TP_TN(D3, hybridTripletExpected, 0.05 > D3_pvalue, net)
         FP_D3bon, FN_D3bon, TP_D3bon, TN_D3bon = calculateD_FP_FN_TP_TN(D3, hybridTripletExpected, (0.05 * 3 / setsOfTriplets) > D3_pvalue, net)
         =#
 
         FP_D, FN_D, TP_D, TN_D = calculate_FP_FN_TP_TN(0.05 > D_pvalue, hybridTripletExpected)
         FP_Dbon, FN_Dbon, TP_Dbon, TN_Dbon = calculate_FP_FN_TP_TN((0.05 * 3 / setsOfTriplets) > D_pvalue, hybridTripletExpected)
 
         FP_Dp, FN_Dp, TP_Dp, TN_Dp = calculate_FP_FN_TP_TN(0.05 > Dp_pvalue, hybridTripletExpected)
         FP_Dpbon, FN_Dpbon, TP_Dpbon, TN_Dpbon = calculate_FP_FN_TP_TN((0.05 * 3 / setsOfTriplets) > Dp_pvalue, hybridTripletExpected)
 
         FP_D3, FN_D3, TP_D3, TN_D3 = calculate_FP_FN_TP_TN(0.05 > D3_pvalue, hybridTripletExpected)
         FP_D3bon, FN_D3bon, TP_D3bon, TN_D3bon = calculate_FP_FN_TP_TN((0.05 * 3 / setsOfTriplets) > D3_pvalue, hybridTripletExpected)
 
         isCorrectDTopology = outgroup_position
         
         push!(df_results, [P1, Hybrid, P2, HyDe_pval, FP, FN, TP, TN, WC, FP_bon, FN_bon, TP_bon, TN_bon, WC_bon,
                            hybridTripletExpected, proposedHybridCorrect, isGhostHybridBool, 
                            isGhostHybridAtRootBool, isMultiHybrid, 
                            lengthToHybridEdge, tipToHybridRatio, 
                            D, D_z, D_pvalue, Dp, Dp_z, Dp_pvalue, D3, D3_z, D3_pvalue, 
                            FP_D, FN_D, TP_D, TN_D, FP_Dbon, FN_Dbon, TP_Dbon, TN_Dbon, 
                            FP_Dp, FN_Dp, TP_Dp, TN_Dp, FP_Dpbon, FN_Dpbon, TP_Dpbon, TN_Dpbon, 
                            FP_D3, FN_D3, TP_D3, TN_D3, FP_D3bon, FN_D3bon, TP_D3bon, TN_D3bon, 
                            true_gamma, proposed_gamma, isCorrectDTopology])
     end
 
     CSV.write(string(outFile), df_results)
 end
 
 #=
 Given separate D-statistics, calculation of true positives, negatives, etc.
 is dependent on whether it matches the species topology required of the D-statistic, 
 as well as if there is a hybrid between the implicated clades in either direction
 
 for positive D statistics, (excess of ABBA) indicates introgression between taxon2 and taxon3
 for negative D statistics, (excess of BABA/"ABAB" in the case of HyDe) 
 indicates introgression between taxon1 and taxon3 in either direction
 
 TODO: method can also be altered to test for whether introgression is between two implicated groups, for now tests if there is a hybrid
 present within the triple tested, given the correct topology.
 =#
 function calculateD_FP_FN_TP_TN(D_statistic, actualHybrid, isSignificant, net)
     
     FP, FN, TP, TN = 0, 0, 0, 0
 
     #=
     t1isHyb = isTrueHybrid(net, triple[1])
     t2isHyb = isTrueHybrid(net, triple[2])
     t3isHyb = isTrueHybrid(net, triple[3])
 
     =#
 
     #=
     if isSignificant
         if D_statistic > 0 # ABBA > BABA 
             if t2isHyb || t3isHyb # AND: these should have a hybrid edge between them, not one caused by ghost hybridization
                 TP = 1
             else
                 FP = 1 
             end
         elseif D_statistic < 0 # BABA > ABBA
             if t1isHyb || t3isHyb # AND: these should have a hybrid edge between them, not one caused by ghost hybridization
                 TP = 1
             else
                 FP = 1
             end
         end
     else # is not significant
         if t1isHyb || t2isHyb || t3isHyb
             FN = 1
         else
             TN = 1
         end
     end
     =#
 
     TP = isSignificant && actualHybrid
     FP = isSignificant && !actualHybrid
     TN = !isSignificant && !actualHybrid
     FN = !isSignificant && actualHybrid
 
     return FP, FN, TP, TN
 
 end
 
 #=
 =#
 function analyzeTICR_MSC(filepath, filename, netname, savefilepath)
     TICR_MSC_Out = DataFrame(CSV.File(filepath))
     net = loadNet(netname)
     outFile = string(savefilepath, filename, "_analyzed.csv")
 
     sig05 = isSignificant(TICR_MSC_Out[!, :MSC_pVal], 0.05)
     sig_bonferroni = isSignificantBonferroni(TICR_MSC_Out[!, :MSC_pVal], 0.05)
 
     setsOfQuartets = size(TICR_MSC_Out,1)
 
     column_names = [:T1, :T2, :T3, :T4, :TICR_pvalue_individual, :TICR_pvalue, :MSC_pvalue,
                      :FP, :FN, :TP, :TN, :FP_bon, :FN_bon, :TP_bon, :TN_bon, 
                     :hybridQuadExpected, :isGhostHybrid, :isGhostAtRoot, 
                     :isMultiHybrid, :lengthToHybridEdge, :tipToHybridRatio, 
                     :true_gamma]
 
     # add all names to dataframe as empty columns
     df_results = DataFrame(column_names .=> Ref([]))
 
     for row in 1:setsOfQuartets
 
         net = loadNet(netname);
         T1 = string(TICR_MSC_Out[row, :t1]); #"1"
         T2 = string(TICR_MSC_Out[row, :t2]);
         T3 = string(TICR_MSC_Out[row, :t3]);
         T4 = string(TICR_MSC_Out[row, :t4]);
         quartet = [T1,T2,T3,T4];
 
         hybridQuadExpected = isHybrid!(net, quartet)
         isGhostHybridBool = false
         isMultiHybrid = false
         isGhostHybridAtRootBool = false
         lengthToHybridEdge, tipToHybridRatio, true_gamma = -1, -1, -1
 
         for i in  1:net.numHybrids
             if i == 2
                 isMultiHybrid = true
             end
             lengthToHybridEdge = calculateLengthToHybridEdge(net, i)
             tipToHybridRatio = calculateTipToHybridRatio(net, i)
             true_gamma = getTrueMixingParameter(net, i)
             isGhostHybridBool = isGhostHybrid(net, i) # is there a bug here?
             isGhostHybridAtRootBool = isGhostHybridAtRoot(net) # this means that there is a hybrid present not ""between"" tested taxa
         end
 
         FP, FN, TP, TN = calculate_FP_FN_TP_TN(sig05[row], hybridQuadExpected)
         FP_bon, FN_bon, TP_bon, TN_bon = calculate_FP_FN_TP_TN(sig_bonferroni[row], hybridQuadExpected)
 
         TICR_pvalue_individual = TICR_MSC_Out[row, :pValTicr]
         TICR_pvalue = TICR_MSC_Out[row, :overallPval]
         MSC_pvalue = TICR_MSC_Out[row, :MSC_pVal]
 
         push!(df_results, [T1, T2, T3, T4, TICR_pvalue_individual, TICR_pvalue, MSC_pvalue, FP, FN, TP, TN,
                      FP_bon, FN_bon, TP_bon, TN_bon, 
                     hybridQuadExpected, isGhostHybridBool, isGhostHybridAtRootBool, 
                     isMultiHybrid, lengthToHybridEdge, tipToHybridRatio, 
                     true_gamma])
     end
 
     CSV.write(string(outFile), df_results)
 end
 
 function main(pathToFile::String, splitby, pathToSaveFile)
    """
     #filename is set up with the format string(netname,"/-gt", num_gene_trees, "-1-1.tre...", HyDe_Dstat.csv
     filename = split(pathToFile, "/")[size(split(pathToFile, "/"), 1)]
     method, netname, netsize, gt_number, trial_number = extractMethodDetails(filename)
     if contains(filename, "-out.txt") && !contains(filename, ".csv") && !contains(filename, "time") && !contains(filename, "50h10")
         analyzeHyDe_DStat(pathToFile, filename, netname, pathToSaveFile)
     elseif contains(filename, "TICR_MSC") && contains(filename, ".csv") && !contains(filename, "time") && !contains(filename, "50h10")
         analyzeTICR_MSC(pathToFile, filename, netname, pathToSaveFile)
     else
         print("input file ", filename, " was not analyzed: does not fit naming conventions")
     end
     """
 end
 

 
 #=
 Tally true positives and negatives from analyzed files creating a file in the form of :
 method     network    networksize    trial   true positives   true negatives   false positives  false negatives   wrong clades
 =#
 function summaryTableTICR(list_output_files, savePath)
     column_names = [:network_name, :net_size, :seq_length, :trial_num, 
     :FP, :FN, :TP, :TN, :FP_bon, :FN_bon, :TP_bon, :TN_bon, :TICR_pval]
     
     summary_table = DataFrame(column_names .=> Ref([]))
 
     for file in list_output_files
         filename = split(file[1], "/")[size(split(file[1], "/"), 1)]
        # print(filename)
         if contains(filename, "analyzed")
             method, netname, netsize, gt_number, trial_number = extractMethodDetails(filename)
             if contains(filename, "TICR")
                 file = DataFrame(CSV.File(file[1]))
 
                 FP = sum(file[!, :FP])
                 FN = sum(file[!, :FN])
                 TP = sum(file[!, :TP])
                 TN = sum(file[!, :TN])
             
                 FP_bon = sum(file[!, :FP_bon])
                 FN_bon = sum(file[!, :FN_bon])
                 TP_bon = sum(file[!, :TP_bon])
                 TN_bon = sum(file[!, :TN_bon])
 
                 TICR_pval = file[1, :TICR_pvalue]
 
                 push!(summary_table, [netname, netsize, gt_number, trial_number,
                     FP, FN, TP, TN, FP_bon, FN_bon, TP_bon, TN_bon, TICR_pval])
             end
         end
     end
     CSV.write(string(savePath), summary_table)
 
 end
 
 
 function extractMethodDetails(filename::AbstractString)
     netname=split(filename, "-")[1]
 
     regex_netsize = r"n\d*"
     netsize = match(regex_netsize, netname)
     netsize = split(netsize.match, "n")[2]
 
     gt_number=split((split(filename, "-gt")[2]), "-")[1]
     regex_trial = r"\d*.tre"
     trial_number = match(regex_trial, filename)
     trial_number = split(trial_number.match, ".tre")[1]
     
     possible_methods = ["HyDe", "D3", "TICR", "MSC"]
     method=nothing
     for met in possible_methods
         if contains(filename, met)
             method = met
         end
     end
 
     return method, netname, netsize, gt_number, trial_number
 end
 
 function summaryTableHyDe(list_output_files, savePath) 
 
     column_names = [:network_name, :net_size, :seq_length, :trial_num, 
                :FP, :FN, :TP, :TN, :WC, :FP_bon, :FN_bon, :TP_bon, :TN_bon, :WC_bon,
                :FP_D, :FN_D, :TP_D, :TN_D, :FP_Dbon, :FN_Dbon, :TP_Dbon, :TN_Dbon, 
                :FP_Dp, :FN_Dp, :TP_Dp, :TN_Dp, :FP_Dpbon, :FN_Dpbon, :TP_Dpbon, :TN_Dpbon, 
                :FP_D3, :FN_D3, :TP_D3, :TN_D3, :FP_D3bon, :FN_D3bon, :TP_D3bon, :TN_D3bon]
     
     summary_table = DataFrame(column_names .=> Ref([]))
 
 
     for file in list_output_files
         filename = split(file[1], "/")[size(split(file[1], "/"), 1)]
         if contains(filename, "analyzed")
             method, netname, netsize, gt_number, trial_number = extractMethodDetails(filename)
             if contains(filename, "-out")
             file = DataFrame(CSV.File(file[1]))
     ### HyDe   
     HyDe_FP = sum(file[!, :FP])
     HyDe_FN = sum(file[!, :FN])
     HyDe_TP = sum(file[!, :TP])
     HyDe_TN = sum(file[!, :TN])
     HyDe_WC = sum(file[!, :WC])
 
     HyDe_FPbon = sum(file[!, :FP_bon])
     HyDe_FNbon = sum(file[!, :FN_bon])
     HyDe_TPbon = sum(file[!, :TP_bon])
     HyDe_TNbon = sum(file[!, :TN_bon])
     HyDe_WCbon = sum(file[!, :WC_bon])
     ### D-statistic
 
    # FP_D,FN_D,TP_D,TN_D,FP_Dbon,FN_Dbon,TP_Dbon,TN_Dbon,FP_Dp,FN_Dp,TP_Dp,TN_Dp,FP_Dpbon,FN_Dpbon,TP_Dpbon,TN_Dpbon,FP_D3,FN_D3,TP_D3,TN_D3,FP_D3bon,FN_D3bon,TP_D3bon,TN_D3bon
     
     Dstat_FP = sum(file[!, :FP_D])
     Dstat_FN = sum(file[!, :FN_D])
     Dstat_TP = sum(file[!, :TP_D])
     Dstat_TN = sum(file[!, :TN_D])
 
     Dstat_FPbon = sum(file[!, :FP_Dbon])
     Dstat_FNbon = sum(file[!, :FN_Dbon])
     Dstat_TPbon = sum(file[!, :TP_Dbon])
     Dstat_TNbon = sum(file[!, :TN_Dbon])
     ### Dp
     Dp_FP = sum(file[!, :FP_Dp])
     Dp_FN = sum(file[!, :FN_Dp])
     Dp_TP = sum(file[!, :TP_Dp])
     Dp_TN = sum(file[!, :TN_Dp])
 
     Dp_FPbon = sum(file[!, :FP_Dpbon])
     Dp_FNbon = sum(file[!, :FN_Dpbon])
     Dp_TPbon = sum(file[!, :TP_Dpbon])
     Dp_TNbon = sum(file[!, :TN_Dpbon])
     ### D3
     D3_FP = sum(file[!, :FP_D3])
     D3_FN = sum(file[!, :FN_D3])
     D3_TP = sum(file[!, :TP_D3])
     D3_TN = sum(file[!, :TN_D3])
 
     D3_FPbon = sum(file[!, :FP_D3bon])
     D3_FNbon = sum(file[!, :FN_D3bon])
     D3_TPbon = sum(file[!, :TP_D3bon])
     D3_TNbon = sum(file[!, :TN_D3bon])
 
     push!(summary_table, [netname, netsize, gt_number, trial_number,
         HyDe_FP, HyDe_FN, HyDe_TP, HyDe_TN, HyDe_WC,
         HyDe_FPbon, HyDe_FNbon, HyDe_TPbon, HyDe_TNbon, HyDe_WCbon,
         Dstat_FP, Dstat_FN, Dstat_TP, Dstat_TN, Dstat_FPbon, Dstat_FNbon, Dstat_TPbon, Dstat_TNbon,
         Dp_FP, Dp_FN, Dp_TP, Dp_TN, Dp_FPbon, Dp_FNbon, Dp_TPbon, Dp_TNbon,
         D3_FP, D3_FN, D3_TP, D3_TN, D3_FPbon, D3_FNbon, D3_TPbon, D3_TNbon])
         end
     end
 end
     
     CSV.write(string(savePath), summary_table)
 end 
 
 function summaryTableHyDe(list_output_files, savePath, filterCondition, filterValue) 
 
     column_names = [:network_name, :net_size, :seq_length, :trial_num, 
                :FP, :FN, :TP, :TN, :WC, :FP_bon, :FN_bon, :TP_bon, :TN_bon, :WC_bon,
                :FP_D, :FN_D, :TP_D, :TN_D, :FP_Dbon, :FN_Dbon, :TP_Dbon, :TN_Dbon, 
                :FP_Dp, :FN_Dp, :TP_Dp, :TN_Dp, :FP_Dpbon, :FN_Dpbon, :TP_Dpbon, :TN_Dpbon, 
                :FP_D3, :FN_D3, :TP_D3, :TN_D3, :FP_D3bon, :FN_D3bon, :TP_D3bon, :TN_D3bon]
     
     summary_table = DataFrame(column_names .=> Ref([]))
 
 
     for file in list_output_files
         filename = split(file[1], "/")[size(split(file[1], "/"), 1)]
         if contains(filename, "analyzed")
             method, netname, netsize, gt_number, trial_number = extractMethodDetails(filename)
             if contains(filename, "HyDe")
             file = DataFrame(CSV.File(file[1]))
     ### HyDe   
     HyDe_FP = sum(file[!, :FP])
     HyDe_FN = sum(file[!, :FN])
     HyDe_TP = sum(file[!, :TP])
     HyDe_TN = sum(file[!, :TN])
     HyDe_WC = sum(file[!, :WC])
 
     HyDe_FPbon = sum(file[!, :FP_bon])
     HyDe_FNbon = sum(file[!, :FN_bon])
     HyDe_TPbon = sum(file[!, :TP_bon])
     HyDe_TNbon = sum(file[!, :TN_bon])
     HyDe_WCbon = sum(file[!, :WC_bon])
     ### D-statistic
 
    # FP_D,FN_D,TP_D,TN_D,FP_Dbon,FN_Dbon,TP_Dbon,TN_Dbon,FP_Dp,FN_Dp,TP_Dp,TN_Dp,FP_Dpbon,FN_Dpbon,TP_Dpbon,TN_Dpbon,FP_D3,FN_D3,TP_D3,TN_D3,FP_D3bon,FN_D3bon,TP_D3bon,TN_D3bon
     
     Dstat_FP = sum(file[!, :FP_D])
     Dstat_FN = sum(file[!, :FN_D])
     Dstat_TP = sum(file[!, :TP_D])
     Dstat_TN = sum(file[!, :TN_D])
 
     Dstat_FPbon = sum(file[!, :FP_Dbon])
     Dstat_FNbon = sum(file[!, :FN_Dbon])
     Dstat_TPbon = sum(file[!, :TP_Dbon])
     Dstat_TNbon = sum(file[!, :TN_Dbon])
     ### Dp
     Dp_FP = sum(file[!, :FP_Dp])
     Dp_FN = sum(file[!, :FN_Dp])
     Dp_TP = sum(file[!, :TP_Dp])
     Dp_TN = sum(file[!, :TN_Dp])
 
     Dp_FPbon = sum(file[!, :FP_Dpbon])
     Dp_FNbon = sum(file[!, :FN_Dpbon])
     Dp_TPbon = sum(file[!, :TP_Dpbon])
     Dp_TNbon = sum(file[!, :TN_Dpbon])
     ### D3
     D3_FP = sum(file[!, :FP_D3])
     D3_FN = sum(file[!, :FN_D3])
     D3_TP = sum(file[!, :TP_D3])
     D3_TN = sum(file[!, :TN_D3])
 
     D3_FPbon = sum(file[!, :FP_D3bon])
     D3_FNbon = sum(file[!, :FN_D3bon])
     D3_TPbon = sum(file[!, :TP_D3bon])
     D3_TNbon = sum(file[!, :TN_D3bon])
 
     push!(summary_table, [netname, netsize, gt_number, trial_number,
         HyDe_FP, HyDe_FN, HyDe_TP, HyDe_TN, HyDe_WC,
         HyDe_FPbon, HyDe_FNbon, HyDe_TPbon, HyDe_TNbon, HyDe_WCbon,
         Dstat_FP, Dstat_FN, Dstat_TP, Dstat_TN, Dstat_FPbon, Dstat_FNbon, Dstat_TPbon, Dstat_TNbon,
         Dp_FP, Dp_FN, Dp_TP, Dp_TN, Dp_FPbon, Dp_FNbon, Dp_TPbon, Dp_TNbon,
         D3_FP, D3_FN, D3_TP, D3_TN, D3_FPbon, D3_FNbon, D3_TPbon, D3_TNbon])
         end
     end
 end
     
     CSV.write(string(savePath), summary_table)
 end 
 
 # julia analyzeMethodOutputFile.jl 
 
 
 
 
 
 
 #=
 For the analysis of HyDe files that don't include a ground truth, 
 but only a network/specified tree to compare to 
 
 
 filepath = file of hyde output to analyze
 filename = desired filename for output to save
 netname = filepath to proposed species tree or network
 savefilepath = desired path to save output
 
 outgroup is unused
 =#
 
 function analyzeHyDe_DStat(filepath, filename, netname, savefilepath, outgroup) 
     HyDeOut = DataFrame(CSV.File(filepath))
     net = loadNet(netname)
     outFile = string(savefilepath, filename, "_analyzed.csv")
 
     sig05 = isSignificant(HyDeOut[!, :Pvalue], 0.05)
     sig_bonferroni = isSignificantBonferroni(HyDeOut[!, :Pvalue], 0.05)
     outgroup = split(filename,"h")[1]
     outgroup = split(outgroup,"n")[1]
     setsOfTriplets = size(HyDeOut,1)
 
     column_names = [:P1, :Hybrid, :P2, :HyDe_pvalue,  
                     :D, :D_z, :D_pvalue, :Dp, :Dp_z, :Dp_pvalue, :D3, :D3_z, :D3_pvalue, 
                     :proposed_gamma, :isDtestTopology]
 
     # add all names to dataframe as empty columns
     df_results = DataFrame(column_names .=> Ref([]))
 
     for row in 1:setsOfTriplets
 
         net = readTopology(netname)
         P1 = string(HyDeOut[row, :P1]);
         Hybrid = string(HyDeOut[row, :Hybrid]);
         P2 = string(HyDeOut[row, :P2]);
         proposed_gamma = HyDeOut[row, :Gamma]
         triplet = [P1,Hybrid,P2];
         quad_includingoutgroup = [P1,Hybrid,P2,outgroup]
 
         ABBA_site = round(HyDeOut[row,:ABBA])
         ABAB_site = round(HyDeOut[row,:ABAB])
         BBAA_site = round(HyDeOut[row,:AABB])
 
         AAAB = round(HyDeOut[row,:AAAB])
         AABA = round(HyDeOut[row,:AABA])
         AABB = round(HyDeOut[row,:AABB])
         AABC = round(HyDeOut[row,:AABC])	
         ABAA = round(HyDeOut[row,:ABAA])	
         ABAB = round(HyDeOut[row,:ABAB])	
         ABAC = round(HyDeOut[row,:ABAC])	
         ABBA = round(HyDeOut[row,:ABBA])	
         BAAA = round(HyDeOut[row,:BAAA])	
         ABBC = round(HyDeOut[row,:ABBC])	
         CABC = round(HyDeOut[row,:CABC])	
         BACA = round(HyDeOut[row,:BACA])	
         BCAA = round(HyDeOut[row,:BCAA])	
         ABCD = round(HyDeOut[row,:ABCD])
 
         # D-statistics should only be calculated where the correct topology is followed. 
         isCorrectDTopology_pos3 = isSisterSisterOutgroup(P1, Hybrid, P2, net)
         isCorrectDTopology_pos1 = isSisterSisterOutgroup(P2, Hybrid, P1, net)
 
         outgroup_position = 0
         if isCorrectDTopology_pos1 
             outgroup_position = 1
         elseif isCorrectDTopology_pos3
             outgroup_position = 3
             # triplet is in the correct order
         else
             outgroup_position = 2
         end
 
         D, D_z, D_pvalue = calcDsd_concatenated(ABBA_site, ABAB_site, BBAA_site, outgroup_position)
         Dp, Dp_z, Dp_pvalue = calcDpsd_concatenated(ABBA_site, ABAB_site, BBAA_site, outgroup_position)
         D3, D3_z, D3_pvalue = calcD3fromHyDe(AAAB, AABA,	AABB,	AABC,	ABAA,	ABAB,	ABAC,	ABBA,	BAAA,	ABBC,	CABC,	BACA,	BCAA,	ABCD, outgroup_position)
         HyDe_pval = HyDeOut[row,:Pvalue]
   
         isCorrectDTopology = outgroup_position
         
         push!(df_results, [P1, Hybrid, P2, HyDe_pval,  
                            D, D_z, D_pvalue, Dp, Dp_z, Dp_pvalue, D3, D3_z, D3_pvalue, 
                            proposed_gamma, isCorrectDTopology])
     end
 
     CSV.write(string(outFile), df_results)
 end
 
 #=
 Remove taxon from trees
 - input tree, taxon to remove, ("Dufourea_novaeangliae")
 - output tree
 =#
 function removeAll(filename, leafname)
 
     newfilename = string(filename, "_remove_", leafname, ".txt")
     io = open(newfilename, "w")
 
     f = readMultiTopology(filename)
 
     for tree in f
         treeString = writeTopology(tree)
         if contains(treeString, leafname)
             PhyloNetworks.deleteleaf!(tree, leafname, keeporiginalroot=true, simplify=true)
         end
         newTopologyString = writeTopology(tree)
         write(io, string(newTopologyString,"\n"))
     end
 
 end