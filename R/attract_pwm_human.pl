#!/usr/bin/perl
# filter a motif PWM file
#Author: Sid

use Getopt::Long;
use Cwd ;

my ($path,$help,$out,$in) = "" ;
GetOptions(
    'help'   => \$help,
    'o=s' => \$out ,
    'i=s' => \$in ,
    'h=s' => \$human_motif,
    'p=s' => \$path ,

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n ";

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	ARGUMENTS :
                    REQUIRED
                    -i -> input attract pwm file
                    -h -> input attract_db metadata file
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }

open(out, "> $out");

# my $human_motif = "ATtRACT_db_human.txt" ;
# my $human_motif = "ATtRACT_db_human_UTR.txt" ;
open (hm, $human_motif)or die "Cannot open file $human_motif \n";
my %hash ;
foreach(<hm>){
	$_ =~ s/\n|\r|\s+$//g ;
  if($_ =~ /^Gene_name/){next ;}
	my(@cols) = split("\t",$_) ;
  my $motif_id = $cols[11] ; my $gene_name = $cols[0] ; my $gene_id = $cols[1] ;
	$hash{$motif_id} = "$gene_name:$gene_id:$motif_id" ;
  #last ;
}
close hm ;

$/ = "\n>" ;
open (in, $in)or die "Cannot open file $in \n";
foreach(<in>) {
	$_ =~ s/\r|\s+$//g ;
  my ($id,$matrix) = split("\n",$_, 2) ;
  $id =~ s/>//;
  $matrix =~ s/\n>//;
  my($motif_id, $motif_len) = split("\t", $id) ;

	if(exists $hash{$motif_id}) {
	 	print ">$hash{$motif_id}\t$motif_id\n$matrix\n" ;
	 }
}

close in ;