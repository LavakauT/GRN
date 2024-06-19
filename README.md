# Informative GRN with kmer (Following KmerDiscovery, ML_pipeline, and kmer_similarity)
><img src="https://github.com/LavakauT/GRN/assets/132649549/9b67a914-730a-4a49-8f11-8fbdd9fdb742" width="40%">


## Prepare files before constructing GRN
1. [VST-transfromed court data](file/vst_norm_count.txt) and Differential group genes ([UU_Gene](file/UU.txt) and [DD_Gene](file/DD.txt) groups in sample data sets).
   You can generate an [original attribute file](attr.txt) before adding Kmer information.
   
2. Distinct and enriched kmer list ([UU_Kmer](training/tss_5utr_promoter_pcc/UU/RF/UU_distinct_pcc_enriched_kmer.txt) and [DD_Kmer](training/tss_5utr_promoter_pcc/DD/RF/DD_distinct_pcc_enriched_kmer.txt)).
  Please follow KmerDiscovery with Differential group genes as TP and non-differential genes ([NN_Gene](file/NN.txt)) as TN to get kmers for motif similarity.

3. Motif similarity table ([UU_Kmer_Similarity](motif/UU/UU_top_sim_dap.txt), [DD_Kmer_Similarity](motif/DD/DD_top_sim_dap.txt)).
   Please follow Kmer_similarity to generate the most similar motif between each kmer and Arabidopsis DAP-seq dataset.

4. [TF_list](Ath_TF_list.txt) and [Gene function table](TAIR10_functional_descriptions.txt). TF list can assign regulators in GRN, and the gene function table can support more information to each node.

## Steps of Constructing GRN
Once you have prepared all the files, please follow the below steps to get the GRN:

1. Follow [ara_network.R](ara_network.R) to generate [GRN_GRNIE3](network/ara_hs_genie3_tf_high.txt) inferred by GENIE3 algorithm and export matrix file/TF list for ARACNe-AP algorithm.
2. Follow [ARACNe-AP](aracne.txt) to install Stand-alone JAVA/ANTs and generate [GRN_ARACNe-AP](network/aracne_network.txt) inferred by ARACNe-AP.
3. Install Cytoscape and upload your two GRN:
```
###### IMPORT GRNs AND GENERATE INTERSECTION GRN ######

> Click Import Network From File System (beside Save Session icon)

> Import GRN_GRNIE3 (network/ara_hs_genie3_tf_high.txt)

> Assign regulatoryGene as Source Node; targetGene as Target Node; Weight as Edge Attribute. Then choose OK.

> Click Import Network From File System (beside the Save Session icon) again.

> Import GRN_ARACNe-AP (network/aracne_network.txt).

> Assign Regulator as Source Node; Target as Target Node; MI/p-value as Edge Attribute. Then choose OK.

> Choose Tools (in tool bar) >> Merge >> Networks. It will jump out Advanced Network Merge window.

> Move two GRNs from the Available network box to the Network to merge box. Choose Intersection and Merge.

> Intersection GRN is done!

###########################
```

```
###### ADD NODE INFORMATION, SET TOPOLOGY FILTER AND EXPORT ######

> In the right panel under the network visualization one, Choose Node Table >> Import Table From Files (besides magnifier icon)

> Import TF_list (Ath_TF_list.txt). Assign Gene_ID as key. Then choose OK.

> Choose Node Table >> Import Table From Files (besides magnifier icon) again.

> Import Gene function table (TAIR10_functional_descriptions.txt). Assign Model_name as key. Then choose OK.

# OPTIONAL: focus on linked nodes
> Choose Filter (The second icon in Layout Tools bar) >> Add new condition >> Topology filter >> Nodes with "less than" >> "3" neighbours within distances "3" >> Apply.
>> Choose Select (in toolbar) >> Edges >> Edges Between Selected Nodes
>>> Choose File (in toolbar) >> New Network >> From Selected Edges, Selected Nodes

> Remember to export Merged Network default edge.csv (network/Merged Network default edge.csv) and Merged Network default node.csv (network/Merged Network default node.csv) from the merged GRN.

> Remember to export Merged GRN file (network/network.cys).

###########################
```

4. Back to [ara_network.R](ara_network.R) and import Merged Network default edge.csv to generate [edge ranking file](network/filter.txt).
5. Back to Cytoscape and load Merged GRN file (network/network.cys).
```
###### ADD EDGE FILTER AND EXPORT SUB-GRN ######

!! TOP20% GRN here is in the row of TOP708 from filter.txt 

> Choose Filter (The second icon in the Layout Tools bar) >> Add new condition >> Column filter >> Edge: MI "is" between "top20 MI filter value" and "highest value" inclusive.

> Choose Filter (The second icon in the Layout Tools bar) >> Add new condition >> Column filter >> Edge: weight "is" between "top20 weight filter value" and "highest value" inclusive.

>> Choose Select (in tool bar) >> Nodes >> Nodes Connected by Selected Edges >> SUB-GRN done!

> Export edge and node file again.

> Remember to save Merged GRN file (network/network.cys).

###########################
```

6. Follow [kmer_application.R](kmer_application.R) and generate [Merged Network default edge_kmer.csv](network/Merged Network default edge_kmer.csv)
7. Back to Cytoscape, load Merged GRN file (network/network.cys), and import Merged Network default edge_kmer.csv in Edge Table.

```
###### ADD KMER INFORMATION IN EDGE AND EXPORT SUB-GRN w/ KMER ######

> In the right panel under the network visualization one, Choose Edge Table >> Import Table From Files (besides magnifier icon)

> Import Merged Network default edge_kmer.csv (network/Merged Network default edge_kmer.csv). Assign name as key. Then choose OK. SUB-GRN w/ kmer done!

> Remember to save Merged GRN file (network/network.cys).

###########################
```

8. Gene centre sub-GRN w/ kmer
```
###### SUBSET SUB-GRN w/ KMER FROM SELECTED TARGET GENE ######

> Choose Filter (The second icon in the Layout Tools bar) >> Add new condition >> Column filter >> Node: name is "your target gene ID"

> Choose Chain (in the bottom, beside Filter) >> choose + >> Add Node Adjacency Transformer >> Take nodes and  "add" "adjacent nodes and edges"

> Choose + below this filter. It will jump out where the "adjacent edges" match the filter

> Choose Column... >> Edge: kmer "is" "Support" (not including Nodata) or Edge: kmer "is not" "Reject" (including Nodata). Repeat and repeat.

> You can export node and edge as new files.

> Remember to save Merged GRN file (network/network.cys).

###########################
```
<img width="461" alt="adjacent edges" src="https://github.com/LavakauT/GRN/assets/132649549/27b2c087-505d-4135-9f60-2e16f250a258">
