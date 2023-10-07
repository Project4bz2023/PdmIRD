#!/usr/bin/perl -w
use strict;
my $dir="/path/to/data";
my $subRVIS="$dir/RVIS.txt"; 
my $DNE="$dir/DNE.txt";
my $GI="$dir/GI.txt";
my $Agscore="$dir/Aggar_score";
my $LOF="$dir/LOF.txt";
my $GVIR="$dir/GVIR";
my $HI="$dir/HI_Predictions";
my $vep_mart="$dir/mart_export.txt.gz";
my $vep_mart1="$dir/mart_export.new";

my $annovar=shift @ARGV;
my $name=`echo -n  $annovar|awk -F"/" '{print \$NF}' |awk -F"." 'BEGIN{OFS=".";ORS=""}{print \$1,\$2,\$3,\$4,\$5}'`;
open SUBRVIS,"$subRVIS" or die $!;
open DNE,"$DNE" or die $!;
open GI,"$GI" or  die $!;
open AG,"$Agscore" or die $!;
open LOF,"$LOF" or die $!;
open GVIR,"$GVIR" or die $!;
open HI,"$HI" or die $!;
open VEP,"gunzip -dc $vep_mart|" or die $!;
open VEP1,"$vep_mart1" or die $!;
open OUT, ">/path/to/result/${name}Ninput" or die $!;
my (%subrvis,%dne,%gi,%ag,%lof,%gvir,%hi,%vep,%vep1);
my ($h_subrvis,$h_dne,$h_gi,$h_ag,$h_lof,$h_gvir,$h_hi);
$h_subrvis=<SUBRVIS>;
chomp($h_subrvis);
$h_subrvis=`echo $h_subrvis|awk 'BEGIN{ORS="";OFS="\t"}{print \$2,\$3}'`;
while(<SUBRVIS>){
	chomp;
	next if (/gene/);
	my ($gene,$other)=(split/\t/,$_,2);
	$subrvis{$gene}=$other;	
}
close SUBRVIS;
while(<VEP>){
        chomp;
        next if (/^gene/i);
        my ($ensg,$gene)=(split/\t/,$_)[0,-1];
	$vep{$ensg}=$gene;
}
close VEP;
while(<VEP1>){
        chomp;
        next if (/^Gene/);
        my ($enst,$nm)=(split/\t/,$_)[2,4];
	if(length($nm) >0){ 
       		$vep1{$nm}=$enst;
	}
}
close VEP1;
$h_dne=<DNE>;
chomp($h_dne);
$h_dne=`echo $h_dne|awk 'BEGIN{ORS="";OFS="\t"}{print \$3,\$4,\$5,\$6,\$7,\$8}'`;
while(<DNE>){
	chomp;
	next if (/gene/);
        my @info=(split/\t/,$_)[0,2..7];
	my $other=join "\t",@info[1..6];
        $dne{$info[0]}=$other;
}
close DNE;
$h_gi=<GI>;
chomp($h_gi);
$h_gi=`echo -n  $h_gi|awk 'BEGIN{ORS="";OFS="\t"}{print \$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$20,\$21,\$22,\$23,\$24}'`;
while(<GI>){
        chomp;
        next if (/gene/i);
        my @info=(split/\t/,$_)[0..16,19..23];
	my $other=join "\t",@info[1..21];
        $gi{$info[0]}=$other;
}
close GI;
$h_ag=<AG>;
chomp($h_ag);
$h_ag=`echo -n  $h_ag|awk 'BEGIN{ORS="";OFS="\t"}{print \$3,\$4}'`;
while(<AG>){
        chomp;
        next if (/gene/);
        my @info=(split/\t/,$_)[0,2,3];
	my $other=join "\t",@info[1,2];
        $ag{$vep{$info[0]}}=$other;
}
close AG;
$h_lof=<LOF>;
chomp($h_lof);
$h_lof=`echo  $h_lof|awk 'BEGIN{ORS="";OFS="\t"}{print \$21,\$22,\$23,\$32,\$33,\$34,\$70}'`;
while(<LOF>){
        chomp;
        next if (/gene/);
	my @info=(split/\t/,$_)[0,20..22,31..33,69];
	my $other=join "\t",@info[1..7];
        $lof{$info[0]}=$other;
}
close LOF;
$h_gvir=<GVIR>;
chomp($h_gvir);
$h_gvir=`echo $h_gvir|awk 'BEGIN{ORS="";OFS="\t"}{print \$2,\$3, \$4,\$5, \$6,\$7}'`;
while(<GVIR>){
	chomp;
	next if (/^gene/);
	my @info=(split/\t/,$_,2);
	$gvir{$info[0]}=$info[1];
}
close GVIR;
$h_hi=<HI>;
chomp($h_hi);
$h_hi=`echo $h_hi|awk 'BEGIN{ORS="";OFS="\t"}{print \$2}'`;
while(<HI>){
        chomp;
        next if (/^gene/);
        my @info=(split/\t/,$_,2);
        $hi{$info[0]}=$info[1];
}
close HI;


open ANNOVAR,"$annovar" or die $!;
my $h_anno=<ANNOVAR>; chomp($h_anno);

print OUT "$h_anno\t$h_subrvis\t$h_dne\t$h_gi\t$h_ag\t$h_lof\t$h_gvir\t$h_hi\n";
while(<ANNOVAR>){
	chomp;
	my ($Tr,$gene) =(split/\t/,$_)[9,6];
	my @AA=(split/;/,$Tr);
	my (@pp,@nm);
	foreach my $i (0..(scalar(@AA)-1)){
		my  @tmp=(split/\:/,$AA[$i])[1,4];
		push @nm,$tmp[0];
		push @pp,$tmp[1];

	}
	if (exists $subrvis{$gene}){
		print OUT "$_\t$subrvis{$gene}\t";
	}
	else{ print OUT "$_\t0.669944\t50.0031\t";  }  ###print mean if the data is not provided
	if (exists $dne{$gene}){
                print OUT "$dne{$gene}\t";
        }
        else{ print OUT "-4.71797\t-5.2589\t-4.89927\t-6.21084\t-6.13192\t-6.0838\t";  }
	if (exists $gi{$gene}){
                print OUT "$gi{$gene}\t";
        }
        else{ print OUT "0.0001\t50.15\t0.0002\t50.08\t0.0001\t50.04\t0.0008\t50.13\t0.00\t50.17\t-0.0001\t50.28\t0.0002\t50.13\t0.0001\t50.25\t50.01\t50.01\t50.01\t5.0e+01\t0.3558\t";  }	
	if (exists $ag{$gene}){
                print OUT "$ag{$gene}\t";
        }
        else{ print OUT "7.93292e-13\t0.00506294\t";  }
	if (exists $lof{$gene}){
                print OUT "$lof{$gene}\t";
        }
        else{ print OUT "0.2451\t0.2568\t0.4981\t-0.20764\t0.75772\t2.1544\t0.3032\t";  }
	if (exists $gvir{$gene}){
                print OUT "$gvir{$gene}\t";
        }
        else{ print OUT "50.0026\t50.0026\t0.991554\t0.983596\t1.00096\t1.00146\t";  }
	if (exists $hi{$gene}){
                print OUT "$hi{$gene}\t";
        }
        else{ print OUT "0.220506\n";  }
}
close ANNOVAR;

