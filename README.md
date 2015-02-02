This program is designed to be a wrapper for analyzing multiple bams with CNV-seq (http://tiger.dbs.nus.edu.sg/cnv-seq/).

cnvSeqWrapper.pl runs cnv-seq commands to analyse all possible pairwise combinations of given bam files.

getCnvSeqRegionsInCommon.pl can then mine the calls for each sample to find regions called in n analyses. 

__INSTALL/RUN__

You need to have samtools installed (preferably version 1.0 or higher). Download the CNV-seq package (http://tiger.dbs.nus.edu.sg/cnv-seq/) and install the R library if you have not done so already. Place the cnv-seq.pl script in the same directory as these scripts and you should be good to go.  

Run cnvSeqWrapper.pl without arguments for usage information. It is assumed you are using a human genome.
 It will generate a folder for each sample containing output files and a folder of [sampleA_vs_sampleB] for each comparison containing an output file with the extension '.cnv.out' containing all the CNV calls for that pairwise comparison. These .cnv.out files can be used to find calls in common between multiple pairwise comparisons using the getCnvSeqRegionsInCommon.pl script. 

Invoke getCnvSeqRegionsInCommon.pl with either --help or --manual flags for help message or manual page. 
