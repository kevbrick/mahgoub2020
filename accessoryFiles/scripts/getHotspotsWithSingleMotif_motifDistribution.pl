#!/usr/bin/perl
#perl getHotspotsWithSingleMotif_motifDistribution.pl --hs [hotpsot BED File] --m [motif PWM file] --g [genome] --dir --out [output File Name]

use strict;
use Getopt::Long;
use Math::Round;

GetOptions (	'hs=s'            => \(my $hs),
				'm=s'             => \(my $motifPWM),
				'g=s'             => \(my $genome),
				'dir+'            => \(my $useDir),
				'out=s'           => \(my $useOut),
				'keepFinalFimo+'  => \(my $keepFinalFimo),
				'initFIMO=s'      => \(my $fimoFolder),
				'plim=s'          => \(my $pLimit = 1),
				'w=i'             => \(my $nWin = 250));

my $randStem  	= "/tmp/tmp_KB_".int(rand()*100000000000);
my $tmpFA  		= $randStem.".fa";
my $tmpErr 		= $randStem.".err";

$fimoFolder 	= $fimoFolder?$fimoFolder:$randStem."_fimo";
my $fimoFile 	= $fimoFolder."/fimo.txt";
my $finalFimo 	= $fimoFolder."/finalfimo.txt";

my $graphFile;
if ($useOut){
		$graphFile = $useOut.'.uniqueHSByPVal.tab';
}else{
		$graphFile = $hs; 
		$graphFile =~ s/^.\///; 
		$graphFile =~ s/^(\S+?)\.\S+$/$1\.uniqueHSByPVal.tab/;
		
}

die if ($graphFile eq $hs);
die ("Please specify the location of perl scripts in the PERLPATH environment variable\n") unless ($ENV{"PERLPATH"});
die ("Please specify the location of perl scripts in the PERLPATH environment variable\n") unless ($ENV{"PERLPATH"});

my $rmList      = "$tmpFA $tmpErr $fimoFile $fimoFolder ";

system("perl ".$ENV{PERLPATH}."/seqFromBed.pl --genome $genome --bed $hs --flanks $nWin --fimo --out $tmpFA >$tmpErr 2>$tmpErr");
system("fimo --max-stored-scores 1000000 --output-pthresh ".($pLimit<5e-3?$pLimit:'5e-3')." --o $fimoFolder $motifPWM $tmpFA") unless (-e $fimoFile);

getOptimalPvalForSingleHSMotif($fimoFile,$finalFimo,$graphFile,$pLimit);

## Keep final motif FIMO file
if ($keepFinalFimo){
	system("cp $finalFimo .");
	
	open FINALFIMOTXT, $finalFimo; 
	open FINALFIMOBED, '>', './finalFimo.tmp.bed'; 
	while (<FINALFIMOTXT>){
		chomp;
		my @F = split(/\t/,$_); 
		my @X = split(/_/,$F[1]); 
		if ($F[2] < $F[3]){
			print FINALFIMOBED join("\t",$X[0],$X[1]+$F[2],$X[1]+$F[3],$F[4].":".$F[5],$F[6].":".$F[7],"+");
		}else{
			print FINALFIMOBED join("\t",$X[0],$X[1]+$F[3],$X[1]+$F[2],$F[4].":".$F[5],$F[6].":".$F[7],"-");
		}
	}
	
	close FINALFIMOBED ; 
	
	system("sort -k1,1 -k2n,2n ./finalFimo.tmp.bed |uniq >finalFimoMotif.bed");
	system("rm ./finalFimo.tmp.bed");
}

my (%hs,%hsV,%hsSc);
open(IN, $finalFimo);

while (<IN>){
	next if ($_ =~ /seque/);
	my @F = split(/\t/,$_);
	my @X = split(/_/,$F[1]); 
	my $strand = ($F[3]>$F[2])?'+':'-';
	
	my ($m);
	if ($useDir){
		#KB Jan 16: $m = ($F[3]>$F[2])?$F[2]:$F[3];
		$m = $F[2];
	}else{
		$m = int(($F[2]+$F[3])/2); 
	}
	
	$hs{$F[1]}++;
	$hsV{$F[1]}  = join("\t",$X[0],($X[1]+$m)-500,($X[1]+$m)+500,$F[4],$F[4],$strand) if ($strand eq '+');
	$hsV{$F[1]}  = join("\t",$X[0],($X[1]+$m)-502,($X[1]+$m)+498,$F[4],$F[4],$strand) if ($strand eq '-');
	$hsSc{$F[1]} = $F[4];
}

my $tf 	= $randStem.".tmpFile";
$rmList .= $tf." ";

open TMP, '>', $tf;

for my $h(keys(%hs)){
	if ($hs{$h} == 1 && $hsV{$h}){
		print TMP $hsV{$h}."\n";
	}
}

if ($useOut){
	system("sort -k1,1 -k2n,2n $tf >$useOut.tab");
}else{
	system("sort -k1,1 -k2n,2n $tf");
}	

####################################################################################
sub getOptimalPvalForSingleHSMotif{
	my ($pfimoFile,$pfinalFimo,$pGraphFile,$pLimz) = @_;
	
	my @testPVals = reverse(5e-3,1e-4,2e-4,4e-4,6e-4,8e-4,1e-5,2e-5,4e-5,6e-5,8e-5,1e-6,1e-7,1e-8);
	
	for my $pv(0..$#testPVals){
		shift @testPVals if ($testPVals[$pv] > $pLimz);
	}
	
	my (%hsCnt);
	
	## Parse FIMO file and count motifs at hotspots for each p-value threshold
	open(IN, $pfimoFile);
	while (<IN>){
		next if ($_ =~ /seque/);
		my ($cs,$hsName,$from,$to,$score,$P,$Q,$seqMatch) = split(/\t/,$_);
		my ($hsCS,$hsFrom,$hsTo) 					 	  = split(/_/,$hsName);
		for my $testP(@testPVals){
			if ($P <= $testP) {
				$hsCnt{$hsName}->{$testP}++;
			}
		}	
	}
	close IN;
	
	## Count number of unique hotspots with a single motif for each P
	my (%UniqueCountsByP,$totalHS);
	for my $hsName(sort keys(%hsCnt)){
		$totalHS++;
		for my $testP(@testPVals){
			next unless ($hsCnt{$hsName}->{$testP});
			$UniqueCountsByP{$testP}++ if ($hsCnt{$hsName}->{$testP} == 1);
		}
	}
	
	## Get best P-value
	my $nMax = -1;
	my $bestP;
	
	open(GRAPH, '>', $pGraphFile);
	
	for my $testP(sort {$a <=> $b} keys(%UniqueCountsByP)){
		print GRAPH join("\t",$testP,$UniqueCountsByP{$testP},sprintf("%4.2f",$UniqueCountsByP{$testP}/$totalHS*100)."\n");
		if ($UniqueCountsByP{$testP} > $nMax){
			$bestP = $testP;
			$nMax  = $UniqueCountsByP{$testP};
		}
	}
	close GRAPH;
	
	plotMe($pGraphFile);
	
	## Write new FIMO
	open(IN, $pfimoFile);
	open(OUT, '>', $pfinalFimo);
	
	while (<IN>){
		next if ($_ =~ /seque/);
		my ($cs,$hsName,$from,$to,$score,$P,$Q,$seqMatch) = split(/\t/,$_);
		my ($hsCS,$hsFrom,$hsTo) 					 	  = split(/_/,$hsName);
		print OUT $_ if ($P <= $bestP && $hsCnt{$hsName}->{$bestP} == 1);
	}
	close IN; close OUT; 
	
}

####################################################################################
sub plotMe{
	my $PMgraphFile = shift;
	
	my $pngFile = $PMgraphFile; $pngFile =~ s/^(.+).tab$/$1.png/;
	my $Rscript = $PMgraphFile; $Rscript =~ s/^(.+).tab$/$1.R/;
	
	#my $Rscript = '/tmp/plotScript_'.int(rand()*10000000000).'.R';
	
	open RS, '>', $Rscript;
	print RS 'library(ggplot2)'."\n";
	print RS 'library(grid)'."\n";
	print RS 'png("'.$pngFile.'",height = 5.5, width = 4.25, units = "in", res=300)'."\n";
	print RS 'myData = read.table("'.$PMgraphFile.'",header=FALSE)'."\n";
	print RS 'colnames(myData)<-c("FIMO P-value","Hotspots (#)","Hotspots (%)")'."\n";
	print RS 'p<-ggplot(myData,aes(x=`FIMO P-value`,y=`Hotspots (%)`))'."\n";
	print RS 'p + geom_line() + scale_x_log10() + geom_point(aes(size=16,guide=FALSE)) + theme(axis.title.x = element_blank())'."\n";
	print RS 'pushViewport(viewport(layout = grid.layout(2, 1)))'."\n";
	print RS 'q<-ggplot(myData,aes(x=`FIMO P-value`,y=`Hotspots (#)`))'."\n";
	print RS 'p1<-p + geom_line() + scale_x_log10() + geom_point(aes(size=16,guide=FALSE)) + theme(axis.title.x = element_blank(), legend.position="none")'."\n";
	print RS 'p2<-q + geom_line(aes(color="red")) + scale_x_log10() + geom_point(aes(size=16,guide=FALSE,color="red")) + theme(legend.position="none")'."\n";
	print RS 'print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))'."\n";
	print RS 'print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))'."\n";
	print RS 'dev.off()'."\n";
	
	close RS; 
	
	system('R --vanilla <'.$Rscript.' >/tmp/R.o 2>/tmp/R.e');
}
