#!/usr/bin/perl
#written by David A. Parry, University of Leeds
#Jan 2015
#
use strict;
use warnings;
use threads;
use File::Basename;
use File::Copy;
use Getopt::Long;
use Cwd;
my @input = ();
my %opts = (input => \@input);
GetOptions(
    \%opts,
    "input=s{,}",    #bam files
    "outdir=s",      
    "pvalue=f",       
    "log2=f",     
    "help",
)or usage("Syntax error");
usage("At least two bams must be provided to the --input option") if @input < 2;
usage() if $opts{help};

sub usage{
    
    my ($message) = @_;
    print "\n$message\n" if $message;

    print <<"EOT"
 
    Usage:
    
    $0 -i <sample1.bam sample2.bam etc.> [options]
    
    Options:
    
    -i, --input  <one or more bam files to analyse>
    -o, --outdir <directory for output files. Defaults to current working directory.>
    -p, --pvalue <p-value to use as a cut-off for cnv-seq.pl. Default is 0.001.>
    -l, --log2   <log2 cut-off for cnv-seq.pl. Default is 0.6.> 
    -h, --help   <this help message>

    This program takes several bam files (at least 2 required) and performs all possible
    pairwise comparisons using cnv-seq.pl (http://tiger.dbs.nus.edu.sg/cnv-seq/). 

    It is assumed that the cnv-seq.pl from http://tiger.dbs.nus.edu.sg/cnv-seq/ is in
    the same directory as this script. You must also have R on your maching and have 
    installed the 'cnv' R library.

    Sample names are inferred from your input file names such that anything preceding the
    first '_' character (if any) is used as the sample name for your file (e.g. for a file
    named sample1_bwa.bam the sample name is inferred to be 'sample1'). A directory is 
    created for each sample where the '.hits' file will be created alongside each .cnv
    and .counts file for each pairwise comparison. 
    
    Final CNV calls are generated in subfolders named in the format 'sample1_vs_sample2' 
    in a file with the extension .cnv.out. A plots folder containing plots for each 
    chromosome containing a CNV call is also created in each of these subfolders. 
    For automated generation of plots for specific CNV calls see the 
    "getCnvSeqRegionsInCommon.pl" script. 
    
    The "getCnvSeqRegionsInCommon.pl" script also provides a means of identifying CNVs 
    identified in multiple pairwise comparisons to help identify 'true' CNVs. 

    This program was written by David A. Parry, University of Leeds, d.a.parry\@leeds.ac.uk
    
    Copyright 2015  David A. Parry

    This program is free software: you can redistribute it and/or modify it under the terms 
    of the GNU General Public License as published by the Free Software Foundation, either 
    version 3 of the License, or (at your option) any later version. This program is 
    distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
    GNU General Public License for more details. You should have received a copy of the 
    GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

EOT
;
    exit 1 if $message;
    exit 0 ;
}
my $original_dir = getcwd();
my $outdir = getcwd();
if ($opts{outdir}){
    $outdir = $opts{outdir};
    if (not -d $outdir){
        print STDERR "Making output directory: $outdir\n";
        mkdir $outdir or die "Could not make output directory $outdir: $!\n";
    }
}

my $pvalue = 0.001;
if ($opts{pvalue}){
    $pvalue = $opts{pvalue};
}

my $log2 = 0.6;
if ($opts{log2}){
    $log2 = $opts{log2};
}

my @samples = ();
my @convert_commands = ();
my @exts = ( qw ( .sam .bam .cram .SAM .BAM .CRAM) ) ;
foreach my $in (@input){
    my ($file, $dir, $ext) = fileparse($in, @exts); 
    my $sample = (split "_", $file)[0];
    push @samples, $sample; 
    print STDERR "Creating directory for $sample...\n"; 
    mkdir "$outdir/$sample" or die "Could not make directory $outdir/$sample: $!\n";
    push @convert_commands, "samtools view -F 4 $in | perl -lane 'print \"\$F[2]\\t\$F[3]\"' > $outdir/$sample/$sample.hits";
}

for (my $i = 0; $i < @convert_commands; $i++){
    my $thr = threads->create(\&runCommand, \@convert_commands, $i);
}

my $err = 0;
while (threads->list(threads::all)){
    foreach my $joinable (threads->list(threads::joinable)){
        $err += $joinable->join();
    }
    sleep 2; 
}
if ($err){
    die "One or more conversion commands failed, exiting...\n";
}

my $n = 0;
foreach my $s (@samples){
    $n++;
    print STDERR "Processing $s (sample $n of " . scalar(@samples) . ")...\n"; 
    chdir("$outdir/$s") or die "Could not change directory to $outdir/$s: $!\n";
    my @others = grep {$s ne $_} @samples;
    die "Need at least two unique samples\n" if not @others;
    my $m = 0; 
    foreach my $o (@others){
        $m++;
        print STDERR "Running cnv-seq.pl for $s vs $o ($m of " . scalar(@others) . ")...\n";
        my $sample_outdir = $s."_vs_". $o;
        mkdir $sample_outdir or die "Could not make directory $sample_outdir: $!\n";
        mkdir "$sample_outdir/plots" or die "Could not make directory $sample_outdir/plots: $!\n";
        #chdir($outdir) or die "Could not change directory to $outdir: $!\n";
    }
    foreach my $o (@others){
        my $sample_outdir = $s."_vs_". $o;
        my $thr = threads->create(\&doCnvSeq, $s, $o, $sample_outdir);
    }
    while (threads->list(threads::all)){
        foreach my $joinable (threads->list(threads::joinable)){
            $joinable->join();
        }
        sleep 2; 
    }
    chdir($original_dir) or die "Could not move directory: $!"; 
}
print STDERR "Finished processing $n samples.\n";

#########################
sub runCommand{
    my ($commands, $i) = @_;
	print STDERR ">>$commands->[$i]\n";
	system "$commands->[$i]";
    return $?;
}

#########################
sub doCnvSeq{
    my ($s, $o, $outdir) = @_;
    system "perl ../cnv-seq.pl --test $s.hits --ref ../$o/$o.hits --genome human";
    die "ERROR: $?\n" if $?;
    print STDERR "Creating and executing R script for $s vs $o...\n";
    my $rscript = make_r_script($s, $o, $outdir);
    system "R CMD BATCH $rscript"; 
    die "ERROR: $?\n" if $?;
    print STDERR "Making chromosome plots for $s vs $o...\n";
    my $plot_script = make_plot_chroms_script($s, $o, $outdir);
    system "R CMD BATCH $plot_script"; 
    die "ERROR: $?\n" if $?;
}
#########################
sub make_plot_chroms_script{
    my ($test, $ref, $outdir) = @_;
    my %chrom = ();
    open (my $IN, "$outdir/$test.hits-vs-$ref.hits.log2-$log2.pvalue-$pvalue.minw-4.cnv.out") 
        or die "Could not open $test.hits-vs-$ref.hits.log2-$log2.pvalue-$pvalue.minw-4.cnv.out: $!";
    while (<$IN>){
        next if /^cnv/;
        my $chr = (split "\t")[1];
        $chr =~ s/^chr//;
        $chrom{$chr}++;
    }
    my $script = "plotChroms-$test-vs-$ref.R"; 
    open (my $PLOT, ">$script" ) or die "Can't create rscript.R: $!\n";
    print $PLOT "library(cnv)\n";
    print $PLOT "data <- read.delim(\"$test.hits-vs-$ref.hits.log2-$log2.pvalue-$pvalue.minw-4.cnv\")\n";
    foreach my $c (keys %chrom){
        print $PLOT "plot.cnv.chr(data, chromosome=\"$c\", ylim = c(-6, 6))\n";
        print $PLOT "ggsave(\"$outdir/plots/chr$c.pdf\")\n";
    }
    close $PLOT;
    return $script; 
}

#########################
sub make_r_script{
    my ($test, $ref, $outdir) = @_;
    my $script = "rscript-$test-vs-$ref.R";
    open (my $R, ">$script") or die "Can't create rscript.R: $!\n";
    
    print $R <<"EOT"
    
    library(cnv)
    data <- read.delim("$test.hits-vs-$ref.hits.log2-$log2.pvalue-$pvalue.minw-4.cnv")
    cnv.print(data, file="$outdir/$test.hits-vs-$ref.hits.log2-$log2.pvalue-$pvalue.minw-4.cnv.out")
    cnv.summary(data)
    plot.cnv.all(data)
    ggsave("$outdir/plots/all_cnvs.pdf")

EOT
;
    close $R;
    return $script;
}
