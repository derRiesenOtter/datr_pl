[1]    9 6793
 Flagging genes and samples with too many missing values...
  ..step 1
pickSoftThreshold: will use block size 6586.
 pickSoftThreshold: calculating connectivity for given powers...
   ..working on genes 1 through 6586 of 6793
   ..working on genes 6587 through 6793 of 6793
   Power SFT.R.sq    slope truncated.R.sq mean.k. median.k. max.k.
1      1 0.955000  2.67000          0.995    3420      3580   4530
2      2 0.907000  1.06000          0.966    2230      2340   3450
3      3 0.756000  0.50300          0.881    1610      1670   2790
4      4 0.313000  0.19700          0.508    1230      1260   2320
5      5 0.000492  0.00675          0.214     982       991   1980
6      6 0.141000 -0.12800          0.345     803       800   1720
7      7 0.321000 -0.23300          0.505     671       653   1520
8      8 0.441000 -0.31900          0.603     570       549   1360
9      9 0.525000 -0.39100          0.672     490       469   1220
10    10 0.577000 -0.45800          0.711     426       404   1100
11    11 0.613000 -0.51500          0.743     374       350   1000
12    12 0.618000 -0.56500          0.758     331       308    921
13    13 0.632000 -0.61500          0.773     295       272    848
14    14 0.649000 -0.65000          0.795     264       240    784
15    15 0.657000 -0.69400          0.804     238       212    728
16    16 0.669000 -0.72900          0.817     216       191    677
17    17 0.682000 -0.76200          0.830     196       171    632
18    18 0.679000 -0.80000          0.828     179       153    593
19    19 0.681000 -0.82800          0.839     164       139    557
20    20 0.690000 -0.84900          0.853     151       126    525
 Calculating module eigengenes block-wise from all genes
   Flagging genes and samples with too many missing values...
    ..step 1
 ....pre-clustering genes to determine blocks..
   Projective K-means:
   ..k-means clustering..
   ..merging smaller clusters...
Block sizes:
gBlocks
   1    2 
4363 2430 
 ..Working on block 1 .
    TOM calculation: adjacency..
    ..will not use multithreading.
     Fraction of slow calculations: 0.000000
    ..connectivity..
    ..matrix multiplication (system BLAS)..
    ..normalization..
    ..done.
   ..saving TOM for block 1 into file lunghemorrhageTOM-block.1.RData
 ....clustering..
 ....detecting modules..
 ....calculating module eigengenes..
 ....checking kME in modules..
     ..removing 33 genes from module 1 because their KME is too low.
     ..removing 12 genes from module 2 because their KME is too low.
     ..removing 4 genes from module 3 because their KME is too low.
 ..Working on block 2 .
    TOM calculation: adjacency..
    ..will not use multithreading.
     Fraction of slow calculations: 0.000000
    ..connectivity..
    ..matrix multiplication (system BLAS)..
    ..normalization..
    ..done.
   ..saving TOM for block 2 into file lunghemorrhageTOM-block.2.RData
 ....clustering..
 ....detecting modules..
 ....calculating module eigengenes..
 ....checking kME in modules..
     ..removing 18 genes from module 1 because their KME is too low.
     ..removing 14 genes from module 2 because their KME is too low.
 ..merging modules that are too close..
     mergeCloseModules: Merging modules whose distance is less than 0.25
       Calculating new MEs...
[1] "Only ME2 seems to be significant"
[1] "Correlations:"
           [,1]
ME1  0.55937083
ME2  0.92768496
ME3  0.23720659
ME6 -0.27880058
ME4  0.13118221
ME5 -0.46302723
ME0  0.02266529
[1] "p-values:"
            [,1]
ME1 0.1173632618
ME2 0.0003118182
ME3 0.5388552395
ME6 0.4675472686
ME4 0.7365605959
ME5 0.2094203898
ME0 0.9538464058
null device 
          1 
