#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

$|=1;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Import SILVAngs export ".../exports/x---xsu---otus.csv
#Output tables for use in MetaPipe. Collapse on taxonomy versions only.

# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("i:m:f:o:ah", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = SILVAngs export otu file (labelled csv but actually tab-delimited)\n";
        print "-m = Location of metapipe directory\n";
        #print "-a = Set to automatically fill in taxonkit output (as done in metapipe)\n";
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
#Expected format: [0] = sample name; [3] = # sequences; [8] = ncbi taxonomic classification; [9] = silva taxonomic classification

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
        if ($split_line[8] ne "")
            {   $ncbi_tax = $split_line[8];
                $ncbi_tax =~ s/^ncbi\|.+\|.+\|//;
                $ncbi_tax =~ s/[^A-Za-z0-9; ]/_/g;    
            }
        else {$ncbi_tax = "Unknown";}
        my $silva_tax;
        if ($split_line[9] ne "")
            {   $silva_tax = $split_line[9];
                $silva_tax =~ s/^silva\|.+\|.+\|//;
                $silva_tax =~ s/[^A-Za-z0-9; ]/_/g;
            }
        else {$silva_tax = "Unknown";}
        my $keytax = $silva_tax.','.$ncbi_tax;
        $TAXA{$keytax}{$sample}{'count'} += $split_line[3];
        unless (exists $TAXA{$keytax}{'asv'})
            {   $TAXA{$keytax}{'asv'} = "ASV_".$asv_count;
                $asv_count += 1;
            }
    }

my @unique_sampleHeaders = uniq @sample_headers;

#Clean up taxonomy
open(NCBI, ">".$options{o}."/ncbiTaxonomy_TOCLEANUP.txt");
foreach my $i (sort keys %TAXA)
    {   my @tax = split(',', $i);
        my @ncbi_taxa = split(';', $tax[1]);
        print NCBI "<$i>";
        foreach my $j (@ncbi_taxa)
            {   print NCBI "\t$j"
            }
        print NCBI "\n";
    }
close(NCBI);

print "\nEdit the NCBI taxonomy file to match a K/P/C/O/F/G/S architecture with tab delimiter\n";
print "Leave the FIRST column AS IS\n";
print "Rename to have _mod.txt at the end of the file name\n";
print "Once finished...press ENTER\n";
my $hold1 = <STDIN>;

#import and set the NCBI taxonomy
open(IN2, "<".$options{o}."/ncbiTaxonomy_TOCLEANUP_mod.txt") or die "\n\nThere is no ncbiTaxonomy_TOCLEANUP_mod.txt in the $options{o} outdirectory!!\n\n";
my @data2 = <IN2>; close(IN2);
foreach my $line (@data2)
	{	chomp($line);
        $line =~ s/ /\+/g;
        $line =~ s/\t/\*/g;
        $line =~ s/\s//;
        $line =~ s/\+/ /g;
        $line =~ s/\*/\t/g;
        $line =~ s/\"//g;
		my @split_line = split('\t', $line);
        my $original_key = $split_line[0];
        $original_key =~ s/\<//;
        $original_key =~ s/\>//;
        my @ntax;
        foreach my $i (1..$#split_line)
            {   if ($split_line[$i] eq "")
                    {   push(@ntax, "NA");   
                    }
                else
                    {   push(@ntax, $split_line[$i]);
                    }
            }
        if (scalar @ntax == 7)
            {   if (exists $TAXA{$original_key})
                    {   my $newncbitax = join(';', @ntax);
                        $TAXA{$original_key}{'cleaned_ncbi'} = $newncbitax;
                    }
                else
                    {   print "\nError matching taxonomy key string $original_key\n";
                    }
            }
        else
            {   if (scalar @ntax > 7)
                    {   die "\n\nYour NCBI modifie taxonomy has extra columns\n\n";}
                if (scalar @ntax < 7)
                    {   if (exists $TAXA{$original_key})
                            {   my @nntax = @ntax;
                                my $term = $#ntax + 1;
                                foreach my $j ($term..6)
                                    {   push(@nntax, "NA");
                                    }
                                @ntax = @nntax;
                                my $newncbitax = join(';', @ntax);
                                $TAXA{$original_key}{'cleaned_ncbi'} = $newncbitax;
                            }
                        else
                            {   print "\nError matching taxonomy key string $original_key\n";
                            }
                    }
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - M A I N - O U T P U T - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
system("mkdir -p ".$options{o}."/ncbi_taxonomy");
system("mkdir -p ".$options{o}."/silva_taxonomy");

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
my @taxTypes = qw | silva ncbi |;
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -