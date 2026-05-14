# run_sctype produces correct snapshots

    Code
      data@meta.data[c("cells_id", "customclassif")]
    Output
                     cells_id             customclassif
      ATGCCAGAACGACT        0       Immune system cells
      CATGGCCTGTGCAT        1       Immune system cells
      GAACCTGATGAACC        2       Immune system cells
      TGACTGGATTCTCA        3       Immune system cells
      AGTCAGACTGCACA        4       Immune system cells
      TCTGATACACGTGT        5       Immune system cells
      TGGTATCTAAACAG        6       Immune system cells
      GCAGCTCTGTTTCT        7       Immune system cells
      GATATAACACGCAT        8       Immune system cells
      AATGTTGACAGTCA        9       Immune system cells
      AGGTCATGAGTGTC       10       Immune system cells
      AGAGATGATCTCGC       11       Immune system cells
      GGGTAACTCTAGTG       12       Immune system cells
      CATGAGACACGGGA       13       Immune system cells
      TACGCCACTCCGAA       14       Immune system cells
      CTAAACCTGTGCAT       15       Immune system cells
      GTAAGCACTCATTC       16       Immune system cells
      TTGGTACTGAATCC       17       Immune system cells
      CATCATACGGAGCA       18       Immune system cells
      TACATCACGCTAAC       19       Immune system cells
      TTACCATGAATCGC       20       Immune system cells
      ATAGGAGAAACAGA       21       Immune system cells
      GCGCACGACTTTAC       22       Immune system cells
      ACTCGCACGAAAGT       23       Immune system cells
      ATTACCTGCCTTAT       24       Immune system cells
      CCCAACTGCAATCG       25       Immune system cells
      AAATTCGAATCACG       26       Immune system cells
      CCATCCGATTCGCC       27 Pancreatic stellate cells
      TCCACTCTGAGCTT       28 Pancreatic stellate cells
      CATCAGGATGCACA       29 Pancreatic stellate cells
      CTAAACCTCTGACA       30 Pancreatic stellate cells
      GATAGAGAAGGGTG       31 Pancreatic stellate cells
      CTAACGGAACCGAT       32 Pancreatic stellate cells
      AGATATACCCGTAA       33 Pancreatic stellate cells
      TACTCTGAATCGAC       34 Pancreatic stellate cells
      GCGCATCTTGCTCC       35 Pancreatic stellate cells
      GTTGACGATATCGG       36 Pancreatic stellate cells
      ACAGGTACTGGTGT       37 Pancreatic stellate cells
      GGCATATGCTTATC       38 Pancreatic stellate cells
      CATTACACCAACTG       39 Pancreatic stellate cells
      TAGGGACTGAACTC       40 Pancreatic stellate cells
      GCTCCATGAGAAGT       41 Pancreatic stellate cells
      TACAATGATGCTAG       42 Pancreatic stellate cells
      CTTCATGACCGAAT       43 Pancreatic stellate cells
      CTGCCAACAGGAGC       44 Pancreatic stellate cells
      TTGCATTGAGCTAC       45 Pancreatic stellate cells
      AAGCAAGAGCTTAG       46 Pancreatic stellate cells
      CGGCACGAACTCAG       47 Pancreatic stellate cells
      GGTGGAGATTACTC       48 Pancreatic stellate cells
      GGCCGATGTACTCT       49 Pancreatic stellate cells
      CGTAGCCTGTATGC       50 Pancreatic stellate cells
      TGAGCTGAATGCTG       51 Pancreatic stellate cells
      CCTATAACGAGACG       52 Pancreatic stellate cells
      ATAAGTTGGTACGT       53 Pancreatic stellate cells
      AAGCGACTTTGACG       54                Mast cells
      ACCAGTGAATACCG       55                Mast cells
      ATTGCACTTGCTTT       56                Mast cells
      CTAGGTGATGGTTG       57                Mast cells
      GCACTAGACCTTTA       58                Mast cells
      CATGCGCTAGTCAC       59                Mast cells
      TTGAGGACTACGCA       60                Mast cells
      ATACCACTCTAAGC       61                Mast cells
      CATATAGACTAAGC       62                Mast cells
      TTTAGCTGTACTCT       63                Mast cells
      GACATTCTCCACCT       64                Mast cells
      ACGTGATGCCATGA       65                Mast cells
      ATTGTAGATTCCCG       66                Mast cells
      GATAGAGATCACGA       67                Mast cells
      AATGCGTGGACGGA       68                Mast cells
      GCGTAAACACGGTT       69                Mast cells
      ATTCAGCTCATTGG       70                Mast cells
      GGCATATGGGGAGT       71                Mast cells
      ATCATCTGACACCA       72                Mast cells
      GTCATACTTCGCCT       73                Mast cells
      TTACGTACGTTCAG       74                Mast cells
      GAGTTGTGGTAGCT       75                Mast cells
      GACGCTCTCTCTCG       76                Mast cells
      AGTCTTACTTCGGA       77                Mast cells
      GGAACACTTCAGAC       78                Mast cells
      CTTGATTGATCTTC       79                Mast cells

# ScTypeAnnotate produces correct annotations

    Code
      res
    Output
      $key
      [1] "this_is_not_an_uuid"
      
      $name
      [1] "ScType-Pancreas-human"
      
      $rootNode
      [1] TRUE
      
      $type
      [1] "cellSets"
      
      $children
      $children[[1]]
      $children[[1]]$key
      [1] "this_is_not_an_uuid"
      
      $children[[1]]$name
      [1] "Immune system cells"
      
      $children[[1]]$rootNode
      [1] FALSE
      
      $children[[1]]$type
      [1] "cellSets"
      
      $children[[1]]$color
      [1] "color_4"
      
      $children[[1]]$cellIds
       [1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
      [26] 25 26
      
      
      $children[[2]]
      $children[[2]]$key
      [1] "this_is_not_an_uuid"
      
      $children[[2]]$name
      [1] "Pancreatic stellate cells"
      
      $children[[2]]$rootNode
      [1] FALSE
      
      $children[[2]]$type
      [1] "cellSets"
      
      $children[[2]]$color
      [1] "color_7"
      
      $children[[2]]$cellIds
       [1] 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51
      [26] 52 53
      
      
      $children[[3]]
      $children[[3]]$key
      [1] "this_is_not_an_uuid"
      
      $children[[3]]$name
      [1] "Mast cells"
      
      $children[[3]]$rootNode
      [1] FALSE
      
      $children[[3]]$type
      [1] "cellSets"
      
      $children[[3]]$color
      [1] "color_1"
      
      $children[[3]]$cellIds
       [1] 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78
      [26] 79
      
      
      

