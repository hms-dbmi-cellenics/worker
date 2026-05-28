# Expression task returns correct values.

    Code
      res
    Output
      $orderedGeneNames
      [1] "MS4A1" "CD79B"
      
      $stats
      $stats$rawMean
      [1] 0.7890259 1.3545382
      
      $stats$rawStdev
      [1] 1.930257 2.203956
      
      $stats$truncatedMin
      [1] 0 0
      
      $stats$truncatedMax
      [1] 5.698965 5.545942
      
      
      $rawExpression
      $rawExpression$values
       [1] 5.779448 4.660719 5.984036 6.666083 5.313305 5.190573 6.106491 5.694729
       [9] 5.537982 5.621232 3.117873 3.449600 4.968821 4.615121 5.779448 5.349125
      [17] 5.697193 6.378826 5.313305 5.190573 5.012325 5.004940 5.537982 6.130608
      [25] 3.378564 3.729311 4.671146 3.330567 4.289671 4.625072 3.831826 3.891037
      [33] 4.641885 4.472554 2.523157
      
      $rawExpression$index
       [1] 10 11 12 13 14 15 16 17 18 19 24 67  0  9 10 11 12 13 14 15 16 17 18 19 29
      [26] 30 32 50 52 53 54 55 58 63 64
      
      $rawExpression$ptr
      [1]  0 12 35
      
      $rawExpression$size
      [1] 80  2
      
      
      $cellIds
      NULL
      

# Expression task keeps order regardless of the request received.

    Code
      res
    Output
      $orderedGeneNames
      [1] "MS4A1" "CD79B"
      
      $stats
      $stats$rawMean
      [1] 0.7890259 1.3545382
      
      $stats$rawStdev
      [1] 1.930257 2.203956
      
      $stats$truncatedMin
      [1] 0 0
      
      $stats$truncatedMax
      [1] 5.698965 5.545942
      
      
      $rawExpression
      $rawExpression$values
       [1] 5.779448 4.660719 5.984036 6.666083 5.313305 5.190573 6.106491 5.694729
       [9] 5.537982 5.621232 3.117873 3.449600 4.968821 4.615121 5.779448 5.349125
      [17] 5.697193 6.378826 5.313305 5.190573 5.012325 5.004940 5.537982 6.130608
      [25] 3.378564 3.729311 4.671146 3.330567 4.289671 4.625072 3.831826 3.891037
      [33] 4.641885 4.472554 2.523157
      
      $rawExpression$index
       [1] 10 11 12 13 14 15 16 17 18 19 24 67  0  9 10 11 12 13 14 15 16 17 18 19 29
      [26] 30 32 50 52 53 54 55 58 63 64
      
      $rawExpression$ptr
      [1]  0 12 35
      
      $rawExpression$size
      [1] 80  2
      
      
      $cellIds
      NULL
      
    Code
      rev_res
    Output
      $orderedGeneNames
      [1] "MS4A1" "CD79B"
      
      $stats
      $stats$rawMean
      [1] 0.7890259 1.3545382
      
      $stats$rawStdev
      [1] 1.930257 2.203956
      
      $stats$truncatedMin
      [1] 0 0
      
      $stats$truncatedMax
      [1] 5.698965 5.545942
      
      
      $rawExpression
      $rawExpression$values
       [1] 5.779448 4.660719 5.984036 6.666083 5.313305 5.190573 6.106491 5.694729
       [9] 5.537982 5.621232 3.117873 3.449600 4.968821 4.615121 5.779448 5.349125
      [17] 5.697193 6.378826 5.313305 5.190573 5.012325 5.004940 5.537982 6.130608
      [25] 3.378564 3.729311 4.671146 3.330567 4.289671 4.625072 3.831826 3.891037
      [33] 4.641885 4.472554 2.523157
      
      $rawExpression$index
       [1] 10 11 12 13 14 15 16 17 18 19 24 67  0  9 10 11 12 13 14 15 16 17 18 19 29
      [26] 30 32 50 52 53 54 55 58 63 64
      
      $rawExpression$ptr
      [1]  0 12 35
      
      $rawExpression$size
      [1] 80  2
      
      
      $cellIds
      NULL
      

