#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);


# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Import ASVs_counts.tsv and x_asvTaxonomyTable.txt
#Output ASV table merged on identical taxonomy

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("a:t:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-a = DADA2 ASVs_counts.tsv file\n";
        print "-t = ASV taxonomy table from processing script (x_asvTaxonomyTable.txt)\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASV;
my %TAX;
my %NEW_ASV;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#IMPORT ASV count table info
open(IN1, "<$options{a}") or die "\n\nThere is no $options{a} file!!\n\n";
my @asv_dat = <IN1>; close(IN1);
my $sample_header = shift(@asv_dat);
my @file_sample_headers = split('\t', $sample_header);
foreach my $line (@asv_dat)
	{	chomp($line);
		my @data = split('\t', $line);
		my $asv = $data[0]; chomp($asv);
		foreach my $i (1..$#data)
			{	my $number = $data[$i]; chomp($number);
                my $samplehead = $file_sample_headers[$i]; chomp($samplehead);
				$ASV{$asv}{$samplehead} = $number; #Stores ASV count by sample
			}
	}
    
#IMPORT ASV Taxonomy table
open(IN2, "<$options{t}") or die "\n\nThere is no $options{t} file!!\n\n";
my @tax_dat = <IN2>; close(IN2); shift(@tax_dat);
foreach my $line (@tax_dat)
	{	chomp($line);
		my @data = split('\t', $line);
		my $asv = $data[0]; chomp($asv);
        my $taxstring = $data[1].";".$data[2].";".$data[3].";".$data[4].";".$data[5].";".$data[6].";".$data[7];
		$TAX{$taxstring}{'asvs'} .= $asv.";";
	}
    
#Collapse ASV counts on identical taxonomy
foreach my $i (sort keys %TAX)
    {   my $common_asvs = $TAX{$i}{'asvs'};
        my @commons = split(';', $common_asvs);
        my $asv_to_use = $commons[0];
        $NEW_ASV{$asv_to_use}{'taxonomy'} = $i;
        foreach my $j (1..$#file_sample_headers)
            {   my $sample = $file_sample_headers[$j]; chomp($sample);
                $NEW_ASV{$asv_to_use}{$sample} = $ASV{$asv_to_use}{$sample};
                foreach my $k (1..$#commons)
                    {   $NEW_ASV{$asv_to_use}{$sample} += $ASV{$commons[$k]}{$sample};
                    }
            }
    }
    
#Print out new ASV table
print "x";
foreach my $i (1..$#file_sample_headers)
    {   my $sample = $file_sample_headers[$i];
        chomp($sample);
        print "\t$sample";
    }
print "\n";
foreach my $i (sort keys %NEW_ASV)
    {   print "$i";
        foreach my $k (1..$#file_sample_headers)
            {   my $sample = $file_sample_headers[$k];
                chomp($sample);
                print "\t$NEW_ASV{$i}{$sample}";
            }
        print "\n";
    }