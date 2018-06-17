normalise
- count the number of reads in each file. Take the largest number of reads & current file number of reads and scale each insertion site. It does not take the difference in the number of insertion sites into account.

The forward and reverse directions are analysed independantly, in addition to a combined analysis. After normalisation, a minimum read threshold is applied.

As this method is reference free, a sliding window approach is used. Fixed windows are generated which can optionally overlap.

For each direction and condition/control, the essentiality is calculated using the windows, independantly.
A comparison is performed between the conditions and controls, for each direction independantly.
The results are filtered, setting a minimum number of reads in the window, and a minimum log fold change, p value and counts per million. This is to remove low level noise. 
The log fold change is then plotted against the genome, where over expression is positve, and under expression is negative. This is perofmred in each direction and the combined total, independantly.  

span gaps


The regions in the original plot files which do notcorrespond to these identified regions as being significat are masked out. 
