#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

$|=1; #Autoflush so pauses within this script are pushed up to the shell script during run.

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Import SILVAngs export ".../exports/x---xsu---otus.csv
#Output tables for use in MetaPipe. Collapse on taxonomy versions only (Nanopore use).
#Create additional option for merged taxonomy from NCBI Eukaryotes.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:m:f:o:ah", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = SILVAngs export otu file (labelled csv but actually tab-delimited)\n";
        print "-m = Location of metapipe directory\n";
        print "-a = Switch to turn on merge of NCBI Eukaryote taxonomy assignments with SILVA Bacteria/Archaeal assignments\n";
        print "-f = Filter cutoff for zzOther\n";
        print "-o = Output directory\n";
        print "-h = This help message\n\n";
        die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %TAXA;
my $asv_count = 1;
my @sample_headers;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
system("mkdir -p ".$options{o});

#Import SILVAngs export file
#Expected format: [0] = sample name; [3] = # sequences; [6] = sequence data; [8] = ncbi taxonomic classification; [9] = silva taxonomic classification

open(IN, "<$options{i}") or die "\n\nThere is no $options{i} file!!\n\n";
my @data = <IN>; close(IN); shift(@data);
my $heads = shift(@data);
if ($heads !~ m/^sample name/)
    {   die "\n\nSILVA export header location has changed and will break this program. Ask for modification.\n\n";
    }
foreach my $line (@data)
	{	chomp($line);
		my @split_line = split('\t', $line);
        my $sample = $split_line[0];
        chomp($sample);
        $sample =~ s/[^A-Za-z0-9_]/_/g;
        $sample = "MP_".$sample;
        push(@sample_headers, $sample);
        my $ncbi_tax;
        my $ncbi_taxid;
        if ($split_line[8] ne "")
            {   $ncbi_tax = $split_line[8];
                $ncbi_tax =~ m/^ncbi\|.+\|(.+)\|/;
                $ncbi_taxid = $1;
                $ncbi_tax =~ s/^ncbi\|.+\|.+\|//;
                $ncbi_tax =~ s/[^A-Za-z0-9; ]/_/g;    
            }
        else {  $ncbi_tax = "Unknown";
                $ncbi_taxid = "NA";
             }
        my $silva_tax;
        if ($split_line[9] ne "")
            {   $silva_tax = $split_line[9];
                $silva_tax =~ s/^silva\|.+\|.+\|//;
                $silva_tax =~ s/[^A-Za-z0-9; ]/_/g;
            }
        else {$silva_tax = "Unknown";}
        my $keytax = $silva_tax.','.$ncbi_tax;
        my $sequencedata = "BLANK";
        if ($split_line[6] ne "")
            {$sequencedata = $split_line[6];}
        $TAXA{$keytax}{$sample}{'count'} += $split_line[3];
        unless (exists $TAXA{$keytax}{'sequence'} && $TAXA{$keytax}{'sequence'} ne "BLANK")
            {   $TAXA{$keytax}{'sequence'} = $sequencedata;
            }
        unless (exists $TAXA{$keytax}{'asv'})
            {   $TAXA{$keytax}{'asv'} = "ASV_".$asv_count;
                $TAXA{$keytax}{'ncbi_taxid'} = $ncbi_taxid;
                $asv_count += 1;
            }
    }

my @unique_sampleHeaders = uniq @sample_headers;

#Clean up taxonomy
my @taxIDS;
foreach my $i (sort keys %TAXA)
    {   unless ($TAXA{$i}{'ncbi_taxid'} eq "NA")
            {   push(@taxIDS, $TAXA{$i}{'ncbi_taxid'});
            }
    }
my @uniq_taxids = uniq @taxIDS;

open(NCBI, ">".$options{o}."/ncbiTaxonomyIDs.txt");
foreach my $i (@uniq_taxids)
    {   print NCBI "$i\n";
    }
close(NCBI);

system("taxonkit lineage ".$options{o}."/ncbiTaxonomyIDs.txt | awk '".'$2!=""'."' > ".$options{o}."/taxonkit_out.txt");
system("taxonkit reformat ".$options{o}."/taxonkit_out.txt | cut -f1,3 > ".$options{o}."/reformatted_taxonkit_out.txt");
system("perl ".$options{m}."/assets/fillIn_taxonkit.pl -i ".$options{o}."/reformatted_taxonkit_out.txt > ".$options{o}."/reformatted_taxonkit_out.txt_temp");
system("mv ".$options{o}."/reformatted_taxonkit_out.txt ".$options{o}."/reformatted_taxonkit_out_ORIGINAL.txt");
system("cat ".$options{o}."/reformatted_taxonkit_out.txt_temp | sed -E 's/[^A-Za-z0-9;[:blank:]]/_/g' > ".$options{o}."/reformatted_taxonkit_out.txt");
system("rm ".$options{o}."/reformatted_taxonkit_out.txt_temp");

open(NCBITAXIN, "<".$options{o}."/reformatted_taxonkit_out.txt") or die "\n\nThere is no reformatted_taxonkit_out.txt in the $options{o} outdirectory!!\n\n";
my @taxiddat = <NCBITAXIN>; close(NCBITAXIN);
my %TAXONKIT;
foreach my $line (@taxiddat)
    {   chomp($line);
        my @splitline = split('\t', $line);
        my $ncbi_taxid = $splitline[0];
        my $ncbi_taxonkit_format;
        if (exists $splitline[1])
            {   $ncbi_taxonkit_format = $splitline[1];
            }
        else
            {   $ncbi_taxonkit_format = "Unknown";
            }
        $TAXONKIT{$ncbi_taxid}{'taxonomystring'} = $ncbi_taxonkit_format;
    }

foreach my $i (sort keys %TAXA)
    {   if (exists $TAXONKIT{$TAXA{$i}{'ncbi_taxid'}}{'taxonomystring'})
            {   my $newNCBItax = $TAXONKIT{$TAXA{$i}{'ncbi_taxid'}}{'taxonomystring'};
                my @splittax = split(';', $newNCBItax);
                if (scalar @splittax < 7)
                    {   my @newlist = @splittax;
                        my $startpos = $#splittax + 1;
                        foreach my $j ($startpos..6)
                            {   push(@newlist, "NA");
                            }
                        @splittax = @newlist;
                    }
                foreach my $j (0..$#splittax)
                    {   if ($splittax[$j] eq "")
                            {   $splittax[$j] = "NA";
                            }
                    }
                my $newnewNCBItax = join(';', @splittax);
                $TAXA{$i}{'cleaned_ncbi'} = $newnewNCBItax;
            }
        else
            {   $TAXA{$i}{'cleaned_ncbi'} = "Unknown;NA;NA;NA;NA;NA;NA";
            }
    }

foreach my $i (sort keys %TAXA)
    {   unless (exists $TAXA{$i}{'cleaned_ncbi'})
            {   die "\n\nSome of the taxa are missing a cleaned NCBI taxonomy\n\n";
            }
    }


open(TAXCHOICES, ">".$options{o}."/silvangs_taxonomyChoices.txt");

foreach my $i (sort keys %TAXA)
    {   my @tax = split(',', $i);
        my @silva_taxa = split(';', $tax[0]);
        my $newtax = "none";
        if (scalar @silva_taxa <= 6)
            {   if ($silva_taxa[0] eq "Bacteria" || $silva_taxa[0] eq "Archaea")
                {   foreach my $j (0..$#silva_taxa)
                        {   unless ($j == 3)
                                {   if ($silva_taxa[$j] =~ m/ales$/)
                                        {   print "\nCheck taxa w/ order-stye name in incorrect column:\n";
                                            print TAXCHOICES "\nCheck taxa w/ order-stye name in incorrect column:\n";
                                            print "$tax[0]\n";
                                            print TAXCHOICES "$tax[0]\n";
                                            print "Input new taxonomy string K/P/C/O/F/G/S with ';' delimiter:\n";
                                            $newtax = <STDIN>;
                                            chomp($newtax);
                                            print TAXCHOICES "Chose:\n";
                                            print TAXCHOICES "$newtax\n";
                                            
                                        }
                                }
                        }
                    foreach my $j (0..$#silva_taxa)
                        {   unless ($j == 4)
                                {   if ($silva_taxa[$j] =~ m/aceae$/)
                                        {   print "\nCheck taxa w/ family-stye name in incorrect column:\n";
                                            print TAXCHOICES "\nCheck taxa w/ family-stye name in incorrect column:\n";
                                            print "$tax[0]\n";
                                            print TAXCHOICES "$tax[0]\n";
                                            print "Input new taxonomy string K/P/C/O/F/G/S with ';' delimiter:\n";
                                            $newtax = <STDIN>;
                                            chomp($newtax);
                                            print TAXCHOICES "Chose:\n";
                                            print TAXCHOICES "$newtax\n";
                                        }
                                }
                        }
                }
            else
                {   print "\nCheck lineage:\n";
                    print TAXCHOICES "\nCheck lineage:\n";
                    print "$tax[0]\n";
                    print TAXCHOICES "$tax[0]\n";
                    print "Input new taxonomy string K/P/C/O/F/G/S with ';' delimiter:\n";
                    $newtax = <STDIN>;
                    chomp($newtax);
                    print TAXCHOICES "Chose:\n";
                    print TAXCHOICES "$newtax\n";
                }
            }
        else
            {   if ($silva_taxa[0] eq "Bacteria" || $silva_taxa[0] eq "Archaea")
                    {   print "\nCheck taxa w/ larger than expected hierarchy for Bacteria/Archaea:\n";
                        print TAXCHOICES "\nCheck taxa w/ larger than expected hierarchy for Bacteria/Archaea:\n";
                        print "$tax[0]\n";
                        print TAXCHOICES "$tax[0]\n";
                        print "Input new taxonomy string K/P/C/O/F/G/S with ';' delimiter:\n";
                        $newtax = <STDIN>;
                        chomp($newtax);
                        print TAXCHOICES "Chose:\n";
                        print TAXCHOICES "$newtax\n";
                    }
                elsif ($silva_taxa[0] eq "Eukaryota")
                    {   print "\nCheck Eukaryota lineage:\n";
                        print TAXCHOICES "\nCheck Eukaryota lineage:\n";
                        print "$tax[0]\n";
                        print TAXCHOICES "$tax[0]\n";
                        print "Input new taxonomy string K/P/C/O/F/G/S with ';' delimiter:\n";
                        $newtax = <STDIN>;
                        chomp($newtax);
                        print TAXCHOICES "Chose:\n";
                        print TAXCHOICES "$newtax\n";
                    }
                else
                    {   print "\nCheck unknown lineage:\n";
                        print TAXCHOICES "\nCheck unknown lineage:\n";
                        print "$tax[0]\n";
                        print TAXCHOICES "$tax[0]\n";
                        print "Input new taxonomy string K/P/C/O/F/G/S with ';' delimiter:\n";
                        print "If you wish to remove this taxa, simply press enter\n";
                        $newtax = <STDIN>;
                        chomp($newtax);
                        print TAXCHOICES "Chose:\n";
                        print TAXCHOICES "$newtax\n";
                    }
            }
        my $final_silva_taxa;
        if ($newtax eq "none")
            {   $final_silva_taxa = join(';', @silva_taxa);
                $final_silva_taxa =~ s/\;uncultured$//;
            }
        elsif ($newtax eq "")
            {   $final_silva_taxa = "DELETE";
            }
        else
            {   my @splitnew = split(';', $newtax);
                @silva_taxa = @splitnew;
                $final_silva_taxa = join(';', @silva_taxa);
            }
        $TAXA{$i}{'cleaned_silva'} = $final_silva_taxa;
    }
close(TAXCHOICES);

foreach my $i (sort keys %TAXA)
    {   unless (exists $TAXA{$i}{'cleaned_silva'})
            {   die "\n\nSome of the taxa are missing a cleaned SILVA taxonomy\n\n";
            }
    }

foreach my $i (sort keys %TAXA)
    {   if ($TAXA{$i}{'cleaned_silva'} eq "DELETE")
            {   delete $TAXA{$i};}
    }

foreach my $i (sort keys %TAXA)
    {   my $silva = $TAXA{$i}{'cleaned_silva'};
        my @split_silvatax = split(';', $silva);
        if (scalar @split_silvatax < 7)
            {   my @nntax = @split_silvatax;
                my $term = $#split_silvatax + 1;
                foreach my $j ($term..6)
                    {   push(@nntax, "NA");
                    }
                @split_silvatax = @nntax;
            }
        foreach my $j (0..$#split_silvatax)
            {   if ($split_silvatax[$j] eq "")
                    {   $split_silvatax[$j] = "NA";
                    }
            }
        my $joinreassign = join(';', @split_silvatax);
        $TAXA{$i}{'cleaned_silva'} = $joinreassign;
    }

if ($options{a}) #Populate $TAXA{$i}{'cleaned_merged'} filling in Euk assignments from NCBI
    {   system("mkdir -p ".$options{o}."/merged_taxonomy");
        open(EUKMERGE, ">".$options{o}."/merged_taxonomy/merged_NCBI_SILVA_Eukaryotes_info.txt");
        print EUKMERGE "ASV\tOriginalSILVA\tNewMergedAssignmet\n";
        foreach my $i (sort keys %TAXA)
            {   my $ncbi_clean = $TAXA{$i}{'cleaned_ncbi'};
                my $silva_clean = $TAXA{$i}{'cleaned_silva'};
                if ($silva_clean =~ m/Bacteria\;Cyanobacteria\;Cyanobacteriia\;Chloroplast/)
                    {   print EUKMERGE "$TAXA{$i}{'asv'}\t";
                        my @split_ncbi = split(';', $ncbi_clean);
                        if ($split_ncbi[0] eq "Eukaryota")
                            {   $TAXA{$i}{'cleaned_merged'} = $ncbi_clean;
                                print EUKMERGE "$silva_clean\t$ncbi_clean\n";
                            }
                        else
                            {   $TAXA{$i}{'cleaned_merged'} = "Eukaryota;NA;NA;NA;NA;NA;NA";
                                print EUKMERGE "$silva_clean\tEukaryota;NA;NA;NA;NA;NA;NA\n";
                            }
                    }
                elsif ($silva_clean =~ m/Bacteria\;Proteobacteria\;Alphaproteobacteria\;Rickettsiales\;Mitochondria/)
                    {   print EUKMERGE "$TAXA{$i}{'asv'}\t";
                        my @split_ncbi = split(';', $ncbi_clean);
                        if ($split_ncbi[0] eq "Eukaryota")
                            {   $TAXA{$i}{'cleaned_merged'} = $ncbi_clean;
                                print EUKMERGE "$silva_clean\t$ncbi_clean\n";
                            }
                        else
                            {   $TAXA{$i}{'cleaned_merged'} = "Eukaryota;NA;NA;NA;NA;NA;NA";
                                print EUKMERGE "$silva_clean\tEukaryota;NA;NA;NA;NA;NA;NA\n";
                            }
                    }
                else
                    {   $TAXA{$i}{'cleaned_merged'} = $silva_clean;
                    }
            }
        close(EUKMERGE);
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - M A I N - O U T P U T - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
system("mkdir -p ".$options{o}."/ncbi_taxonomy");
system("mkdir -p ".$options{o}."/silva_taxonomy");

my @taxTypes;
if ($options{a})
    {   @taxTypes = qw | silva ncbi merged |;
        
    }
else
    {   @taxTypes = qw | silva ncbi |;
        
    }


#ASVs_counts.tsv = Plain counts for all "ASVs". In this case, each ASV simply represents a unique taxa.
open(ASV_count, ">".$options{o}."/ASVs_counts.tsv");
print ASV_count "x";
foreach my $i (0..$#unique_sampleHeaders)
    {   print ASV_count "\t$unique_sampleHeaders[$i]";
    }
print ASV_count "\n";

foreach my $i (sort keys %TAXA)
    {   my $chosen_asv = $TAXA{$i}{'asv'};
        print ASV_count "$chosen_asv";
        foreach my $j (0..$#unique_sampleHeaders)
            {   if (exists $TAXA{$i}{$unique_sampleHeaders[$j]})
                    {   my $print_value = $TAXA{$i}{$unique_sampleHeaders[$j]}{'count'};
                        print ASV_count "\t$print_value";
                    }
                else
                    {   print ASV_count "\t0";
                    }
            }
        print ASV_count "\n";
    }
close(ASV_count);

#ASVs_counts_NOUNKNOWNS.tsv = Plain counts for all "ASVs", except those with unknown taxonomy.
foreach my $taxonomyType (0..$#taxTypes)
    {   open(ASV_count, ">".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/ASVs_counts_NOUNKNOWNS.tsv");
        print ASV_count "x";
        foreach my $i (0..$#unique_sampleHeaders)
            {   print ASV_count "\t$unique_sampleHeaders[$i]";
            }
        print ASV_count "\n";
        foreach my $i (sort keys %TAXA)
            {   my $currentTaxonomy = $TAXA{$i}{'cleaned_'.$taxTypes[$taxonomyType]};
                unless ($currentTaxonomy eq "Unknown;NA;NA;NA;NA;NA;NA")
                    {   my $chosen_asv = $TAXA{$i}{'asv'};
                        print ASV_count "$chosen_asv";
                        foreach my $j (0..$#unique_sampleHeaders)
                            {   if (exists $TAXA{$i}{$unique_sampleHeaders[$j]})
                                    {   my $print_value = $TAXA{$i}{$unique_sampleHeaders[$j]}{'count'};
                                        print ASV_count "\t$print_value";
                                    }
                                else
                                    {   print ASV_count "\t0";
                                    }
                            }
                        print ASV_count "\n";
                    }
            }
        close(ASV_count);
    }

#morphology_asvTaxonomyTable.txt = ASV to taxonomy.
foreach my $taxonomyType (0..$#taxTypes)
    {   open(ASV_taxa, ">".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/silvangs_".$taxTypes[$taxonomyType]."_asvTaxonomyTable.txt");
        print ASV_taxa "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
        foreach my $i (sort keys %TAXA)
            {   print ASV_taxa $TAXA{$i}{'asv'};
                my $clean_taxa = $TAXA{$i}{'cleaned_'.$taxTypes[$taxonomyType]};
                my @clean_taxa_array = split(';', $clean_taxa);
                foreach my $j (0..$#clean_taxa_array)
                    {   print ASV_taxa "\t".$clean_taxa_array[$j];
                    }
                print ASV_taxa "\n";
            }
        close(ASV_taxa);
    }

#morphology_asvTaxonomyTable_NOUNKNOWNS.txt = ASV to taxonomy (without UNKNOWN).
foreach my $taxonomyType (0..$#taxTypes)
    {   open(ASV_taxa, ">".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/silvangs_".$taxTypes[$taxonomyType]."_asvTaxonomyTable_NOUNKNOWNS.txt");
        print ASV_taxa "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
        foreach my $i (sort keys %TAXA)
            {   my $clean_taxa = $TAXA{$i}{'cleaned_'.$taxTypes[$taxonomyType]};
                unless ($clean_taxa eq "Unknown;NA;NA;NA;NA;NA;NA")
                    {   print ASV_taxa $TAXA{$i}{'asv'};
                        my @clean_taxa_array = split(';', $clean_taxa);
                        foreach my $j (0..$#clean_taxa_array)
                            {   print ASV_taxa "\t".$clean_taxa_array[$j];
                            }
                        print ASV_taxa "\n";
                    }
            }
        close(ASV_taxa);
    }

#KRONA Plots
foreach my $taxonomyType (0..$#taxTypes)
    {   system("mkdir -p ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/KRONA_inputs");
        foreach my $uniquesample (0..$#unique_sampleHeaders)
            {   open(ASV_taxa, ">".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/KRONA_inputs/".$unique_sampleHeaders[$uniquesample]."_KRONA.txt");
                foreach my $i (sort keys %TAXA)
                    {   if (exists $TAXA{$i}{$unique_sampleHeaders[$uniquesample]})
                            {   my $print_value = $TAXA{$i}{$unique_sampleHeaders[$uniquesample]}{'count'};
                                print ASV_taxa "$print_value";
                                my $tax = $TAXA{$i}{'cleaned_'.$taxTypes[$taxonomyType]};
                                my @tax_array = split(';', $tax);
                                foreach my $entry (0..$#tax_array)
                                    {   if ($tax_array[$entry] ne "NA")
                                            {   print ASV_taxa "\t$tax_array[$entry]";
                                            }
                                        else
                                            {   last;
                                            }
                                    }
                                print ASV_taxa "\n";
                            }
                    }
                close(ASV_taxa);
            }
            
            open(ASV_taxa, ">".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/KRONA_inputs/silvangs_samplesSummedKRONA.txt");
            foreach my $uniquesample (0..$#unique_sampleHeaders)
            {   foreach my $i (sort keys %TAXA)
                    {   if (exists $TAXA{$i}{$unique_sampleHeaders[$uniquesample]})
                            {   my $print_value = $TAXA{$i}{$unique_sampleHeaders[$uniquesample]}{'count'};
                                print ASV_taxa "$print_value";
                                my $tax = $TAXA{$i}{'cleaned_'.$taxTypes[$taxonomyType]};
                                my @tax_array = split(';', $tax);
                                foreach my $entry (0..$#tax_array)
                                    {   if ($tax_array[$entry] ne "NA")
                                            {   print ASV_taxa "\t$tax_array[$entry]";
                                            }
                                        else
                                            {   last;
                                            }
                                    }
                                print ASV_taxa "\n";
                            }
                    }
            }
            close(ASV_taxa);
    }

foreach my $taxonomyType (0..$#taxTypes)
    {   my @locationKRONA;
        foreach my $i (@unique_sampleHeaders)
            {   push(@locationKRONA, $options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/KRONA_inputs/".$i."_KRONA.txt");
            }
        my $printkronasamples = join(' ', @locationKRONA);
        system("ImportText.pl -o ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/silvangs_master_krona.html $printkronasamples");
        
        system("ImportText.pl -o ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/silvangs_samplesSummedKRONA.html ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/KRONA_inputs/silvangs_samplesSummedKRONA.txt");
    }

#ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv = relative abundance ASV biom table
foreach my $taxonomyType (0..$#taxTypes)
    {   open(ASV_count, ">".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv");
        print ASV_count "x";
        foreach my $i (0..$#unique_sampleHeaders)
            {   print ASV_count "\t$unique_sampleHeaders[$i]";
            }
        print ASV_count "\n";
        
        my @sample_sums;
        foreach my $j (0..$#unique_sampleHeaders)
            {   my $samp_sum = 0;
                foreach my $i (sort keys %TAXA)
                    {   my $clean_taxa = $TAXA{$i}{'cleaned_'.$taxTypes[$taxonomyType]};
                        unless ($clean_taxa eq "Unknown;NA;NA;NA;NA;NA;NA")
                            {   if (exists $TAXA{$i}{$unique_sampleHeaders[$j]})
                                            {   my $value = $TAXA{$i}{$unique_sampleHeaders[$j]}{'count'};
                                                $samp_sum += $value;
                                            }
                            }
                    }
                push(@sample_sums, $samp_sum);
            }
        
        foreach my $i (sort keys %TAXA)
            {   my $clean_taxa = $TAXA{$i}{'cleaned_'.$taxTypes[$taxonomyType]};
                unless ($clean_taxa eq "Unknown;NA;NA;NA;NA;NA;NA")
                    {   my $chosen_asv = $TAXA{$i}{'asv'};
                        print ASV_count "$chosen_asv";
                        foreach my $j (0..$#unique_sampleHeaders)
                            {   if (exists $TAXA{$i}{$unique_sampleHeaders[$j]})
                                    {   my $value = $TAXA{$i}{$unique_sampleHeaders[$j]}{'count'};
                                        my $relabund = 100 * $value / $sample_sums[$j];
                                        print ASV_count "\t$relabund";
                                    }
                                else
                                    {   print ASV_count "\t0";
                                    }
                            }
                        print ASV_count "\n";
                    }
            }
        close(ASV_count);
    }
    
open(SAMP, ">".$options{o}."/sample_order_fromSILVAfile.txt");
foreach my $i (sort @unique_sampleHeaders)
    {   print SAMP "$i\n";
    }
close(SAMP);

#Run filter_lowabundance_taxa.pl on outfiles to create zzOther taxonomy file
foreach my $taxonomyType (0..$#taxTypes)
    {   system('perl '.$options{m}.'/assets/filter_lowabundance_taxa.pl -a '.$options{o}."/".$taxTypes[$taxonomyType].'_taxonomy/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv -t '.$options{o}."/".$taxTypes[$taxonomyType].'_taxonomy/silvangs_'.$taxTypes[$taxonomyType].'_asvTaxonomyTable_NOUNKNOWNS.txt -p '.$options{f}.' > '.$options{o}.'/'.$taxTypes[$taxonomyType].'_taxonomy/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt');
    }


#Move KRONA plots to Figures directory
foreach my $taxonomyType (0..$#taxTypes)
    {   system("mkdir -p ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/Figures/00_KRONA_plots");
        system("cp ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/silvangs_master_krona.html ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/Figures/00_KRONA_plots/silvangs_".$taxTypes[$taxonomyType]."Taxonomy_master_krona.html");
        system("cp ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/KRONA_plots/silvangs_samplesSummedKRONA.html ".$options{o}."/".$taxTypes[$taxonomyType]."_taxonomy/Figures/00_KRONA_plots/silvangs_".$taxTypes[$taxonomyType]."Taxonomy_samplesSummedKRONA.html");
    }

open(HEADS, ">".$options{o}."/unique_TaxaFolders.txt");
foreach my $i (@taxTypes)
    {   print HEADS "${i}\n";
    }
close(HEADS);

open(SEQDATA, ">".$options{o}."/ASVs.fa");
foreach my $i (sort keys %TAXA)
    {   print SEQDATA ">$TAXA{$i}{'asv'}\n";
        print SEQDATA "$TAXA{$i}{'sequence'}\n";
    }
close(SEQDATA);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -