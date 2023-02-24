# run_sctype produces correct snapshots

    Code
      data@meta.data[c("cells_id", "customclassif")]
    Output
         cells_id             customclassif
      1         0       Immune system cells
      2         1       Immune system cells
      3         2       Immune system cells
      4         3       Immune system cells
      5         4       Immune system cells
      6         5       Immune system cells
      7         6       Immune system cells
      8         7       Immune system cells
      9         8       Immune system cells
      10        9       Immune system cells
      11       10       Immune system cells
      12       11       Immune system cells
      13       12       Immune system cells
      14       13       Immune system cells
      15       14       Immune system cells
      16       15       Immune system cells
      17       16       Immune system cells
      18       17       Immune system cells
      19       18       Immune system cells
      20       19       Immune system cells
      21       20       Immune system cells
      22       21       Immune system cells
      23       22       Immune system cells
      24       23       Immune system cells
      25       24       Immune system cells
      26       25       Immune system cells
      27       26       Immune system cells
      28       27 Pancreatic stellate cells
      29       28 Pancreatic stellate cells
      30       29 Pancreatic stellate cells
      31       30 Pancreatic stellate cells
      32       31 Pancreatic stellate cells
      33       32 Pancreatic stellate cells
      34       33 Pancreatic stellate cells
      35       34 Pancreatic stellate cells
      36       35 Pancreatic stellate cells
      37       36 Pancreatic stellate cells
      38       37 Pancreatic stellate cells
      39       38 Pancreatic stellate cells
      40       39 Pancreatic stellate cells
      41       40 Pancreatic stellate cells
      42       41 Pancreatic stellate cells
      43       42 Pancreatic stellate cells
      44       43 Pancreatic stellate cells
      45       44 Pancreatic stellate cells
      46       45 Pancreatic stellate cells
      47       46 Pancreatic stellate cells
      48       47 Pancreatic stellate cells
      49       48 Pancreatic stellate cells
      50       49 Pancreatic stellate cells
      51       50 Pancreatic stellate cells
      52       51 Pancreatic stellate cells
      53       52 Pancreatic stellate cells
      54       53 Pancreatic stellate cells
      55       54                Mast cells
      56       55                Mast cells
      57       56                Mast cells
      58       57                Mast cells
      59       58                Mast cells
      60       59                Mast cells
      61       60                Mast cells
      62       61                Mast cells
      63       62                Mast cells
      64       63                Mast cells
      65       64                Mast cells
      66       65                Mast cells
      67       66                Mast cells
      68       67                Mast cells
      69       68                Mast cells
      70       69                Mast cells
      71       70                Mast cells
      72       71                Mast cells
      73       72                Mast cells
      74       73                Mast cells
      75       74                Mast cells
      76       75                Mast cells
      77       76                Mast cells
      78       77                Mast cells
      79       78                Mast cells
      80       79                Mast cells

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
      
      
      

