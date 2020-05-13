#!/usr/bin/perl -w
#use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);


# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Take ASV count table, taxAssignX.txt, reformatted_taxonkit_out, genbank common names from NCBI.
#Required settings: filtering options for species and genus cutoffs, output name.
#Optional inputs: list of ASVs to ignore, list of samples in the order desired for output.
#Outputs: Reformatted out with one ASV per sample per line, Krona plots, Bar charts (with and without unknowns),
#         multi-ASV taxa heatmap, common names outputs, list of unknown ASVs. All with and without filtering.


# - - - - - C O M M A N D    L I N E    O P T I O N S - - - - - - - -
my %options=();
getopts("a:s:t:n:f:c:d:o:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-a = ASV counts table (make sure there is text in the upperleft)\n";
        print "-s = DADA2 taxonomy output (w/ headers ASV Perc Len TaxID correction)\n";
	print "-t = reformatted taxonkit output (TaxID\tsimplified taxonomy)\n";
	print "-f = filtering options Species,Genus (e.g. 97,95)\n";
	print "-n = Allin Output basename\n";
	print '-c = Location of common names file (grep "genbank common name" from names.dmp NCBI taxonomy file)'; print "\n";
	print "-d = List of ASVs to ignore (one per line) for outputs ignoring contaminants and/or unknowns\n";
	print "-o = List of samples (one per line) in the order you want them exported. Must be exact matches to ASV counts table.\n";
	print "       Does not have to include all samples.\n";
	print "-h = This help message\n\n";
	die;
    }

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %ASV;
my %TAXON;
my @uniq_bartaxa_ig;
my @sample_headers;
my @file_sample_headers;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#IMPORT ASV count table info
open(IN1, "<$options{a}") or die "\n\nThere is no $options{a} file!!\n\n";
my @asv_dat = <IN1>; close(IN1);
my $sample_header = shift(@asv_dat);
@file_sample_headers = split('\t', $sample_header);
my @new_array_headers;
foreach my $userheaders (@file_sample_headers)
	{	my $value = $userheaders;
		chomp($userheaders);
		push(@new_array_headers, $userheaders);	
	}
@file_sample_headers = @new_array_headers;
foreach my $line (@asv_dat)
	{	chomp($line);
		my @data = split('\t', $line);
		my $asv = $data[0]; chomp($asv);
		foreach my $i (1..$#data)
			{	my $number = $data[$i];
				chomp($number);
				$ASV{$asv}{$file_sample_headers[$i]} = $number; #Stores ASV count by sample
			}
	}

#IMPORT SAMPLE ORDER (if applicable)
if ($options{o})
	{	open(SAMP, "<$options{o}") or die "\n\nThere is no $options{o} file!!\n\n";
		my @usersampleorder = <SAMP>; close(SAMP);
		foreach my $i (@usersampleorder)
			{	chomp($i);
				push(@sample_headers, $i);
			}
	}
else
	{	shift(@file_sample_headers);
		@sample_headers = @file_sample_headers;
	}

foreach my $i (sort keys %ASV)
	{	foreach my $j (0..$#sample_headers)
			{	if (exists $ASV{$i}{$sample_headers[$j]})
					{}
				else {die "\n\nImported sample headers don't match the ASV count table. Try again.\n\n";}
			}
	}


#IMPORT DADA2 taxonomy output
open(IN2, "<$options{s}") or die "\n\nThere is no $options{s} file!!\n\n";
my @dada_dat = <IN2>; close(IN2); shift(@dada_dat);
foreach my $line (@dada_dat)
	{	chomp($line);
		my @data = split('\t', $line);
		my $blah = $data[0]; chomp($blah);
		$ASV{$blah}{'perchit'} = $data[1];
		$ASV{$blah}{'lenhit'} = $data[2];
		$ASV{$blah}{'taxid'} = $data[3];
		$ASV{$blah}{'correction'} = $data[4];
		my $filterinfo = $options{f};
		my @splitfilter = split(',', $filterinfo);
		if ($data[1] >= $splitfilter[0])
			{	$ASV{$blah}{'confidence'} = "SPECIES";	
			}
		if ($data[1] < $splitfilter[0] && $data[1] >= $splitfilter[1])
			{	$ASV{$blah}{'confidence'} = "GENUS";
			}
		if ($data[1] < $splitfilter[1] && $data[1] >= 90)
			{	$ASV{$blah}{'confidence'} = "ORDER";
			}
		if ($data[1] < 90)
			{	$ASV{$blah}{'confidence'} = "NOCONFIDENCE";
			}
	}

#IMPORT reformatted taxonkit taxonomy out
open(IN3, "<$options{t}") or die "\n\nThere is no $options{t} file!!\n\n";
my @taxon_dat = <IN3>; close(IN3);
foreach my $lines (@taxon_dat)
	{	chomp($lines);
		my @taxon_split = split('\t', $lines);
		my $taxonclean = $taxon_split[0];
		chomp($taxonclean);
		my $cleantaxastring = $taxon_split[1];
		chomp($cleantaxastring);
		$TAXON{$taxonclean}{'taxastring'} = $cleantaxastring;
	}

#IMPORT common names
open(INT, "<$options{c}") or die "\n\nThere is no $options{c} file!!\n\n";
my @commondat = <INT>; close(INT);
foreach my $linet (@commondat)
	{	chomp($linet);
		my @common_split = split('\t', $linet);
		my $taxid_clean = $common_split[0];
		chomp($taxid_clean);
		my $common_name_clean = $common_split[2];
		chomp($common_name_clean);
		if (exists $TAXON{$taxid_clean})
			{	$TAXON{$taxid_clean}{'common_name'} = $common_name_clean;
			}
	}
	
#IMPORT List of ASVs to ingore in output
my @ignore_array;
my $ignore_ASV_string;
if ($options{d})
	{	open(IGNORE, "<$options{d}") or die "\n\nThere is no $options{d} file!!\n\n";
		my @igdat = <IGNORE>; close(IGNORE);
		foreach my $linei (@igdat)
			{	chomp($linei);
				push(@ignore_array, $linei);
			}
		foreach my $i (@ignore_array)
			{	$ignore_ASV_string .= "_".$i."_";
			}
	}

	
#Basic filtering of taxonomy to:
#1) replace ;;;;;;;string with Environmental unknown
#2) replace and double ;; with ;@@; for later trimming and detection
#3) replace Eukarya;; with Environmental unknown

foreach my $i (sort keys %TAXON)
	{	my $taxastring = $TAXON{$i}{'taxastring'};
		if ($taxastring =~ m/^;/)
			{	$TAXON{$i}{'taxastring'} = "Environmental Unknown;@@;@@;@@;@@;@@;@@";
				$taxastring = $TAXON{$i}{'taxastring'};
			}
		if ($taxastring =~ m/;;/)
			{	my @splitstring = split(';', $taxastring);
				if ($#splitstring == 6)
					{	foreach my $j (0..$#splitstring)
							{	if ($splitstring[$j] eq '')
									{	$splitstring[$j] = "@@";
										my $next = $j + 1;
										foreach my $k ($next..$#splitstring)
											{	$splitstring[$k] = "@@";
											}
										last;	
									}
							}
						my $newgoodats = join(';', @splitstring);
						if ($newgoodats !~ m/^Eukaryota;@/)
							{	$TAXON{$i}{'taxastring'} = $newgoodats;
							}
						else {$TAXON{$i}{'taxastring'} = "Environmental Unknown;@@;@@;@@;@@;@@;@@";}
					}
				else {die "\nYou forgot to load the reformatted taxonkit file with K/P/C/O/F/G/S only!\n\n";}
			}
	}

#CHOOSE TAXONOMY for output
foreach my $i (sort keys %ASV)
	{	if (exists $ASV{$i}{'taxid'})	#CREATE list of taxids to consider
		{	my $asv_taxid = $ASV{$i}{'taxid'};
			my @tax_list;
			if ($asv_taxid =~ m/\;/ || $asv_taxid =~ m/\,/)
				{	my @multi = split(',', $asv_taxid);
					foreach my $j (@multi)
						{	$j =~ s/\ //;
							if ($j =~ m/\;/)
								{	my @multipli = split(';', $j);
									foreach my $k (@multipli)
										{	push(@tax_list, $k);	
										}
								}
							else {push(@tax_list, $j);}
						}
				}
			else {push(@tax_list, $asv_taxid);}
			my $choiceindicator = 0;
			my $commonnameYES = 0;
			my $desiredcommonname;
			my $choice; # SET the first taxid in the list as the first choice
			if ($ASV{$i}{'confidence'} eq "SPECIES") # SET depth of first choice based on confidence (%)
				{	my $tap = $TAXON{$tax_list[0]}{'taxastring'};
					if ($tap !~ m/@@/)
						{	$choice = $TAXON{$tax_list[0]}{'taxastring'};
							$choiceindicator = 1;
							if (exists $TAXON{$tax_list[0]}{'common_name'})
								{	$commonnameYES = 1;
									$desiredcommonname = $TAXON{$tax_list[0]}{'common_name'};
								}
						}
					else
						{	my $newbie = $TAXON{$tax_list[0]}{'taxastring'};
							my @split_newbie = split(';', $newbie);
							my @pleasework;
							foreach my $position (0..$#split_newbie)
								{	if ($split_newbie[$position] eq "@@")
										{	last;
										}
									else {
										push(@pleasework, $split_newbie[$position]);
									}	
								}
							my $newwithout = join(';', @pleasework);
							$choice = $newwithout;
						}
				}
			if ($ASV{$i}{'confidence'} eq "GENUS")
				{	my $firstchoice = $TAXON{$tax_list[0]}{'taxastring'};
					chomp($firstchoice);
					my @splitmodtax = split(';', $firstchoice);
					pop(@splitmodtax);
					if ($splitmodtax[$#splitmodtax] ne "@@")
						{	my $newmodtax = join(';', @splitmodtax);
							$choice = $newmodtax;
							$choiceindicator = 1;
						}
					else
						{	my @pleasework;
							foreach my $position (0..$#splitmodtax)
								{	if ($splitmodtax[$position] eq "@@")
										{	last;
										}
									else {
										push(@pleasework, $splitmodtax[$position]);
									}	
								}
							my $newwithout = join(';', @pleasework);
							$choice = $newwithout;
						}
				}
			if ($ASV{$i}{'confidence'} eq "ORDER")
				{	my $firstchoice = $TAXON{$tax_list[0]}{'taxastring'};
					chomp($firstchoice);
					my @splitmodtax = split(';', $firstchoice);
					pop(@splitmodtax); pop(@splitmodtax); pop(@splitmodtax);
					if ($splitmodtax[$#splitmodtax] ne "@@")
						{	my $newmodtax = join(';', @splitmodtax);
							$choice = $newmodtax;
							$choiceindicator = 1;
						}
					else
						{	my @pleasework;
							foreach my $position (0..$#splitmodtax)
								{	if ($splitmodtax[$position] eq "@@")
										{	last;
										}
									else {
										push(@pleasework, $splitmodtax[$position]);
									}	
								}
							my $newwithout = join(';', @pleasework);
							$choice = $newwithout;
						}
				}
			if ($ASV{$i}{'confidence'} eq "NOCONFIDENCE")
				{	$choice = "Unknown";
				}
			if ($#tax_list > 0) # SET COMPARISON taxonomy if taxid list has more than one
				{	foreach my $heck (1..$#tax_list)
						{	my $compare;
							my $commonnameCOMP = 0;
							my $compindicator = 0;
							if ($ASV{$i}{'confidence'} eq "SPECIES")
								{	my $tap = $TAXON{$tax_list[$heck]}{'taxastring'};
									if ($tap !~ m/@@/)
										{	$compare = $TAXON{$tax_list[$heck]}{'taxastring'};
											$compindicator = 1;
											if (exists $TAXON{$tax_list[$heck]}{'common_name'})
												{	$commonnameCOMP = 1;
												}
										}
									else
										{	if ($choiceindicator == 0)	# any ;; would end up here
												{	my $newbie = $TAXON{$tax_list[$heck]}{'taxastring'};
													my @split_newbie = split(';', $newbie);
													my @pleasework;
													foreach my $position (0..$#split_newbie)
														{	if ($split_newbie[$position] eq "@@")
																{	last;
																}
															else {
																push(@pleasework, $split_newbie[$position]);
															}	
														}
													my $newwithout = join(';', @pleasework);
													$compare = $newwithout;
												}	
											if ($choiceindicator == 1)	# only if no wonky tax assignment
												{	$compare = $choice;
												}
											
										}
								}
							if ($ASV{$i}{'confidence'} eq "GENUS")
								{	my $firstchoice = $TAXON{$tax_list[$heck]}{'taxastring'};
									chomp($firstchoice);
									my @splitmodtax = split(';', $firstchoice);
									pop(@splitmodtax);
									if ($splitmodtax[$#splitmodtax] ne "@@")
										{	my $newmodtax = join(';', @splitmodtax);
											$compare = $newmodtax;
											$compindicator = 1;
										}
									else
										{	if ($choiceindicator == 0)
												{	my @pleasework;
													foreach my $position (0..$#splitmodtax)
														{	if ($splitmodtax[$position] eq "@@")
																{	last;
																}
															else {
																push(@pleasework, $splitmodtax[$position]);
															}	
														}
													my $newwithout = join(';', @pleasework);
													$compare = $newwithout;
												}
											if ($choiceindicator == 1)
												{	$compare = $choice;
												}
										}
								}
							if ($ASV{$i}{'confidence'} eq "ORDER")
								{	my $firstchoice = $TAXON{$tax_list[$heck]}{'taxastring'};
									chomp($firstchoice);
									my @splitmodtax = split(';', $firstchoice);
									pop(@splitmodtax); pop(@splitmodtax); pop(@splitmodtax);
									if ($splitmodtax[$#splitmodtax] ne "@@")
										{	my $newmodtax = join(';', @splitmodtax);
											$compare = $newmodtax;
											$compindicator = 1;
										}
									else
										{	if ($choiceindicator == 0)
												{	my @pleasework;
													foreach my $position (0..$#splitmodtax)
														{	if ($splitmodtax[$position] eq "@@")
																{	last;
																}
															else {
																push(@pleasework, $splitmodtax[$position]);
															}	
														}
													my $newwithout = join(';', @pleasework);
													$compare = $newwithout;
												}
											if ($choiceindicator == 1)
												{	$compare = $choice;	
												}
										}
								}
							if ($ASV{$i}{'confidence'} eq "NOCONFIDENCE")
								{	$compare = "Unknown";
								}
							my @new_choice;
							#print "$choice\t$compare\n";
							unless ($compindicator == 1 && $choiceindicator == 0)
								{
									if ($choice =~ m/;/ && $compare =~ m/;/)
										{	my @choice_array = split(';', $choice);
											my @compare_array = split(';', $compare);
											foreach my $s (0..$#choice_array)
												{	if (exists $compare_array[$s])
														{	if ($choice_array[$s] eq $compare_array[$s])
																{push(@new_choice, $choice_array[$s]);
																}
															else {last;}
														}
												}
										}
									if ($choice =~ m/;/ && $compare !~ m/;/)
										{	my @choice_array = split(';', $choice);
											my @compare_array;
											push(@compare_array,$compare);
											foreach my $s (0)
												{	if ($choice_array[$s] eq $compare_array[$s])
														{push(@new_choice, $choice_array[$s]);
														}
													else {last;}
												}
										}
									if ($choice !~ m/;/ && $compare =~ m/;/)
										{	my @choice_array;
											push(@choice_array, $choice);
											my @compare_array = split(';', $compare);
											foreach my $s (0)
												{	if ($choice_array[$s] eq $compare_array[$s])
														{push(@new_choice, $choice_array[$s]);
														}
													else {last;}
												}
										}
									if ($choice !~ m/;/ && $compare !~ m/;/)
										{	if ($choice eq $compare)
												{	push(@new_choice, $choice);
												}
										
										}
									if (scalar(@new_choice) == 0)
										{$choice = "Unknown";}
									else {	$choice = join(';', @new_choice);
										if ($#new_choice != 6)
											{	$commonnameYES = 0;
											}
										if ($#new_choice == 6)
											{	if ($commonnameCOMP == 1 && $commonnameYES == 1)
													{	my $choice_common = $desiredcommonname;
														my $compare_common = $TAXON{$tax_list[$heck]}{'common_name'};
														if ($choice_common ne $compare_common)
															{	$commonnameYES = 0;
															}
														
													}
												if ($commonnameCOMP == 1 && $commonnameYES == 0)
													{	$desiredcommonname = $TAXON{$tax_list[$heck]}{'common_name'};
														$commonnameYES = 1;
													}
											}
									      }
								}
							if ($compindicator == 1 && $choiceindicator == 0)
								{	$choice = $compare;
									$choiceindicator = 1;
									if ($commonnameCOMP == 1)
										{	$desiredcommonname = $TAXON{$tax_list[$heck]}{'common_name'};
											$commonnameYES = 1;
										}
								}
							
						}
				}
			$ASV{$i}{'finaltaxachoice'} = $choice;
			if ($commonnameYES == 1)
				{	$ASV{$i}{'common_name'} = $desiredcommonname;
				}
		}
		else {$ASV{$i}{'finaltaxachoice'} = "Unknown";}
	}


##Outputs
foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		open($sample, ">".$sample."_KRONA.txt");
		#print $sample "count\ttaxa\n";
	}

if ($options{d}) {
foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		my $sample_ig = $sample."_IGNORE";
		open($sample_ig, ">".$sample_ig."_KRONA.txt");
		#print $sample "count\ttaxa\n";
	}
}

open(OUT, ">".$options{n}."_allin_KRONA.txt");
open(WHOLEKRONA, ">".$options{n}."_wholeKRONA.txt");
open(ASVTAX, ">".$options{n}."_asvTaxonomyTable.txt");
print ASVTAX "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
print OUT "Sample\tASV\tcount\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t(common name)\n";
foreach my $i (sort keys %ASV)
	{   my $taxput = $ASV{$i}{'finaltaxachoice'};
		my @taxputs = split(';', $taxput);
        print ASVTAX "$i";
        foreach my $l (@taxputs) {
			print ASVTAX "\t$l";}
		if ($#taxputs < 6)
			{	my $actualtaxdepth = $#taxputs + 1;
				foreach my $printtab ($actualtaxdepth..6)
					{	print ASVTAX "\tNA";
					}
			}
        print ASVTAX "\n";
        foreach my $j (0..$#sample_headers)
			{	my $hello = $sample_headers[$j];
				chomp($hello);
				unless ($ASV{$i}{$sample_headers[$j]} == 0)
				{
				print OUT "$sample_headers[$j]\t";
				print OUT "$i\t";
				print OUT "$ASV{$i}{$sample_headers[$j]}";
				foreach my $l (@taxputs) {
				print OUT "\t$l";}
				if ($#taxputs < 6)
					{	my $actualtaxdepth = $#taxputs + 1;
						foreach my $printtab ($actualtaxdepth..6)
							{	print OUT "\t";
							}
					}
				if (exists $ASV{$i}{'common_name'})
					{	print OUT "\t(".$ASV{$i}{'common_name'}.")";
					}
				else {print OUT "\t";}
				print OUT "\n";
				print $hello "$ASV{$i}{$sample_headers[$j]}";
				print WHOLEKRONA "$ASV{$i}{$sample_headers[$j]}";
				foreach my $l (@taxputs) {print $hello "\t$l"; print WHOLEKRONA "\t$l";}
				if (exists $ASV{$i}{'common_name'})
					{	print $hello " (".$ASV{$i}{'common_name'}.")";
						print WHOLEKRONA " (".$ASV{$i}{'common_name'}.")";
					}
				print $hello "\n";
				print WHOLEKRONA "\n";
				}
			}
	}
close(OUT);
close(WHOLEKRONA);
close(ASVTAX);

foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		close($sample);
	}

if ($options{d}) {
open(OUT_IG, ">".$options{n}."_IGNORE_allin_KRONA.txt");
print OUT_IG "Sample\tASV\tcount\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t(common name)\n";
open(ASVTAX_IG, ">".$options{n}."_IGNORE_asvTaxonomyTable.txt");
print ASVTAX_IG "ASV\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
foreach my $i (sort keys %ASV)
	{	unless ($ignore_ASV_string =~ m/_${i}_/)
		{
        
        my $taxput = $ASV{$i}{'finaltaxachoice'};
		my @taxputs = split(';', $taxput);
        print ASVTAX_IG "$i";
        foreach my $l (@taxputs) {
			print ASVTAX_IG "\t$l";}
		if ($#taxputs < 6)
			{	my $actualtaxdepth = $#taxputs + 1;
				foreach my $printtab ($actualtaxdepth..6)
					{	print ASVTAX_IG "\t";
					}
			}
        print ASVTAX_IG "\n";
		foreach my $j (0..$#sample_headers)
			{	my $hello = $sample_headers[$j];
				chomp($hello);
				my $hello_ig = $hello."_IGNORE";
				unless ($ASV{$i}{$sample_headers[$j]} == 0)
				{
				print OUT_IG "$sample_headers[$j]\t";
				print OUT_IG "$i\t";
				print OUT_IG "$ASV{$i}{$sample_headers[$j]}";
				foreach my $l (@taxputs) {
				print OUT_IG "\t$l";}
				if ($#taxputs < 6)
					{	my $actualtaxdepth = $#taxputs + 1;
						foreach my $printtab ($actualtaxdepth..6)
							{	print OUT_IG "\t";
							}
					}
				if (exists $ASV{$i}{'common_name'})
					{	print OUT_IG "\t(".$ASV{$i}{'common_name'}.")";
					}
				else {print OUT_IG "\t";}
				print OUT_IG "\n";
				print $hello_ig "$ASV{$i}{$sample_headers[$j]}";
				foreach my $l (@taxputs) {print $hello_ig "\t$l";}
				if (exists $ASV{$i}{'common_name'})
					{	print $hello_ig " (".$ASV{$i}{'common_name'}.")";
					}
				print $hello_ig "\n";
				}
			}
		}
	}
close(OUT_IG);
close(ASVTAX_IG);

foreach my $i (0..$#sample_headers)
	{	my $sample = $sample_headers[$i];
		chomp($sample);
		my $sample_ig = $sample."_IGNORE";
		close($sample_ig);
	}
}

#BARCHART OUT
my @bartaxa;
my @commonnamebartaxa;
foreach my $j (0..$#sample_headers)
	{	foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0)
					{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
						chomp($pushchoice);
						if ($pushchoice =~ m/;/)
							{	my @splittingthis = split(';', $pushchoice);
								my $lastassignment = $splittingthis[$#splittingthis];
								chomp($lastassignment);
								if (exists $ASV{$i}{'common_name'})
									{	my $newlastassignment = $lastassignment." (".$ASV{$i}{'common_name'}.")";
										push(@commonnamebartaxa, $newlastassignment);
									}
								push(@bartaxa, $lastassignment);
							}
						else {push(@bartaxa, $pushchoice);}
					}
			}
	}
my @uniq_bartaxa = uniq @bartaxa;

if ($options{d})
	{ my @bartaxa_ig;
	  foreach my $j (0..$#sample_headers)
		{	foreach my $i (sort keys %ASV)
				{	unless ($ASV{$i}{$sample_headers[$j]} == 0 || $ignore_ASV_string =~ m/_${i}_/)
						{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
							chomp($pushchoice);
							if ($pushchoice =~ m/;/)
								{	my @splittingthis = split(';', $pushchoice);
									my $lastassignment = $splittingthis[$#splittingthis];
									chomp($lastassignment);
									push(@bartaxa_ig, $lastassignment);
								}
							else {push(@bartaxa_ig, $pushchoice);}
						}
				}
		}
	@uniq_bartaxa_ig = uniq @bartaxa_ig;
	}


open(TAXAOUT, ">".$options{n}."_unique_terminaltaxa.txt");
foreach my $i (@uniq_bartaxa)
	{	unless($i eq "Unknown" || $i eq "Environmental Unknown" || $i =~ m/__/){
		print TAXAOUT "$i\n";}
	}
close(TAXAOUT);
my $namme = $options{n};
system("taxonkit name2taxid ".$namme."_unique_terminaltaxa.txt > ".$namme."_unique_terminaltaxa_w_taxids.txt");

open(TID, "<".$options{n}."_unique_terminaltaxa_w_taxids.txt") or die "\n\nSomething wrong with taxonkit name2taxid output\n\n";
my @commonname_last_dat = <TID>; close(TID);
my %CommonName_Term;
foreach my $line (@commonname_last_dat)
	{	chomp($line);
		my @line_split = split('\t', $line);
		if (exists $line_split[1])
            {my $taxid_clean = $line_split[1];
            chomp($taxid_clean);
            $CommonName_Term{$taxid_clean}{'name'} = $line_split[0];
            }
	}
open(COT, ">".$options{n}."_taxid_to_commonname_ALL.txt");
foreach my $linet (@commondat)
	{	chomp($linet);
		my @common_split = split('\t', $linet);
		my $taxid_clean = $common_split[0];
		chomp($taxid_clean);
		my $common_name_clean = $common_split[2];
		chomp($common_name_clean);
		if (exists $CommonName_Term{$taxid_clean})
			{	print COT "$CommonName_Term{$taxid_clean}{'name'} (".$common_name_clean.")\t$taxid_clean\n";
			}
	}



unless (scalar(@commonnamebartaxa) == 0)
	{	my @uniq_commonnamebartaxa = uniq @commonnamebartaxa;
		open(COMMONBAR, ">".$options{n}."_SPECIESlevel_commonNames_for_barchart.txt");
		foreach my $i (@uniq_commonnamebartaxa)
			{	print COMMONBAR "$i\n";
			}
		close(COMMONBAR);
	}

open(BARCHART, ">".$options{n}."_barchart.txt");
open(BARCHART_forR, ">".$options{n}."_barchart_forR.txt");
print BARCHART_forR "Value\tSample\tTerminalTaxa\n";
print BARCHART "Sample\t";
foreach my $i (@uniq_bartaxa)
	{	print BARCHART "$i\t";
	}
print BARCHART "\n";
foreach my $j (0..$#sample_headers)
	{	my $barsampleheader = $sample_headers[$j];
		chomp($barsampleheader);
		print BARCHART "$barsampleheader\t";
		my %BARCHART;
		foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0)
					{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
						chomp($pushchoice);
						if ($pushchoice =~ m/;/)
							{	my @splittingthis = split(';', $pushchoice);
								my $lastassignment = $splittingthis[$#splittingthis];
								chomp($lastassignment);
								$BARCHART{$lastassignment} += $ASV{$i}{$sample_headers[$j]};
							}
						else {$BARCHART{$pushchoice} += $ASV{$i}{$sample_headers[$j]};}
					}
			}
		foreach my $uniq_taxa (@uniq_bartaxa)
			{	my $hit = 0;
				foreach my $k (sort keys %BARCHART)
					{	if ($uniq_taxa eq $k)
							{	$hit = 1;
								print BARCHART "$BARCHART{$k}\t";
                                print BARCHART_forR "$BARCHART{$k}\t$sample_headers[$j]\t$uniq_taxa\n";
							}
					}
				if ($hit == 0)
					{	print BARCHART "0\t";
					}
			}
		print BARCHART "\n";
	}
close(BARCHART);


if ($options{d})
{
open(BARCHART_ig, ">".$options{n}."_IGNORE_barchart.txt");
print BARCHART_ig "Sample\t";
foreach my $i (@uniq_bartaxa_ig)
	{	print BARCHART_ig "$i\t";
	}
print BARCHART_ig "\n";
foreach my $j (0..$#sample_headers)
	{	my $barsampleheader = $sample_headers[$j];
		chomp($barsampleheader);
		print BARCHART_ig "$barsampleheader\t";
		my %BARCHART;
		foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0 || $ignore_ASV_string =~ m/_${i}_/)
					{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
						chomp($pushchoice);
						if ($pushchoice =~ m/;/)
							{	my @splittingthis = split(';', $pushchoice);
								my $lastassignment = $splittingthis[$#splittingthis];
								chomp($lastassignment);
								$BARCHART{$lastassignment} += $ASV{$i}{$sample_headers[$j]};
							}
						else {$BARCHART{$pushchoice} += $ASV{$i}{$sample_headers[$j]};}
					}
			}
		foreach my $uniq_taxa (@uniq_bartaxa_ig)
			{	my $hit = 0;
				foreach my $k (sort keys %BARCHART)
					{	if ($uniq_taxa eq $k)
							{	$hit = 1;
								print BARCHART_ig "$BARCHART{$k}\t";
							}
					}
				if ($hit == 0)
					{	print BARCHART_ig "0\t";
					}
			}
		print BARCHART_ig "\n";
	}
close(BARCHART_ig);
}


#run KRONA plot
my @kronasamples;
foreach my $i (0..$#sample_headers)
	{	my $hebs = $sample_headers[$i];
		chomp($hebs);
		my $newhebs = $hebs."_KRONA.txt";
		push(@kronasamples, $newhebs);
	}
my $printkronasamples = join(' ', @kronasamples);
my $naming = $options{n};
system("ImportText.pl -o ".$naming."_master_krona.html $printkronasamples");
my $wholekronaout = $options{n}."_wholeKRONA.txt";
system("ImportText.pl -o ".$naming."_wholeKRONA.html $wholekronaout");

if ($options{d})
	{	my @kronasamples_ig;
		foreach my $i (0..$#sample_headers)
			{	my $hebs = $sample_headers[$i];
				chomp($hebs);
				my $newhebs = $hebs."_IGNORE_KRONA.txt";
				push(@kronasamples_ig, $newhebs);
			}
		my $printkronasamples_ig = join(' ', @kronasamples_ig);
		my $naming = $options{n};
		system("ImportText.pl -o ".$naming."_IGNORE_master_krona.html $printkronasamples_ig");
	}

##Barchart without unknowns
open(NOUNKNOWN, ">".$options{n}."_NO_UNKNOWNS_barchart.txt");
print NOUNKNOWN "Sample\t";
foreach my $i (@uniq_bartaxa)
	{	unless ($i eq "Unknown" || $i eq "Environmental Unknown")
			{	print NOUNKNOWN "$i\t";
			}
	}
print NOUNKNOWN "\n";
foreach my $j (0..$#sample_headers)
	{	my $barsampleheader = $sample_headers[$j];
		chomp($barsampleheader);
		print NOUNKNOWN "$barsampleheader\t";
		my %BARCHART;
		foreach my $i (sort keys %ASV)
			{	unless ($ASV{$i}{$sample_headers[$j]} == 0)
					{	my $pushchoice = $ASV{$i}{'finaltaxachoice'};
						chomp($pushchoice);
						if ($pushchoice =~ m/;/)
							{	my @splittingthis = split(';', $pushchoice);
								my $lastassignment = $splittingthis[$#splittingthis];
								chomp($lastassignment);
								$BARCHART{$lastassignment} += $ASV{$i}{$sample_headers[$j]};
							}
						else {$BARCHART{$pushchoice} += $ASV{$i}{$sample_headers[$j]};}
					}
			}
		foreach my $uniq_taxa (@uniq_bartaxa)
			{	unless ($uniq_taxa eq "Unknown" || $uniq_taxa eq "Environmental Unknown")
					{	my $hit = 0;
						foreach my $k (sort keys %BARCHART)
							{	if ($uniq_taxa eq $k)
									{	$hit = 1;
										print NOUNKNOWN "$BARCHART{$k}\t";
									}
							}
						if ($hit == 0)
							{	print NOUNKNOWN "0\t";
							}
					}
			}
		print NOUNKNOWN "\n";
	}
close(NOUNKNOWN);



##Taxa with shared ASV w/ heatmap?
my %TAXHEAT;

foreach my $i (sort keys %ASV)
	{	my $tax = $ASV{$i}{'finaltaxachoice'};
		chomp($tax);
		if ($tax =~ m/;/)
			{	my @split_tax = split(';', $tax);
				foreach my $j (0..$#split_tax)
					{	if (exists $TAXHEAT{$split_tax[$j]})
							{	$TAXHEAT{$split_tax[$j]}{'asvs'} .= ";".$i;
							}
						else
							{	$TAXHEAT{$split_tax[$j]}{'asvs'} = $i;
								$TAXHEAT{$split_tax[$j]}{'depth'} = $j + 1;
							}
					}
			}
		else
			{	if (exists $TAXHEAT{$tax})
					{	$TAXHEAT{$tax}{'asvs'} .= ";".$i;
					}
				else
					{	$TAXHEAT{$tax}{'asvs'} = $i;
						$TAXHEAT{$tax}{'depth'} = 1;
					}
			}	
	}

open(OUTDEEP1, ">".$options{n}."_heatmap_multiASV.txt");
foreach my $i (sort keys %TAXHEAT)
	{	if ($TAXHEAT{$i}{'asvs'} =~ m/;/)
			{	print OUTDEEP1 "$i <<Depth ".$TAXHEAT{$i}{'depth'}.">>\tSample\t";
				my $multiasv = $TAXHEAT{$i}{'asvs'};
				my @multiasv_split = split(';', $multiasv);
				foreach my $j (@multiasv_split)
					{	print OUTDEEP1 "$j\t";
					}
				print OUTDEEP1 "\n";
				foreach my $k (0..$#sample_headers)
					{	print OUTDEEP1 "$i <<Depth ".$TAXHEAT{$i}{'depth'}.">>\t";
						print OUTDEEP1 "$sample_headers[$k]\t";
						foreach my $j (@multiasv_split)
							{	print OUTDEEP1 "$ASV{$j}{$sample_headers[$k]}\t";
							}
						print OUTDEEP1 "\n";
					}
			}
	}
close(OUTDEEP1);


open(UNKNOWNS, ">".$options{n}."_unknown_asvids.txt");
foreach my $i (sort keys %ASV)
	{	my $taxonunknown = $ASV{$i}{'finaltaxachoice'};
		chomp($taxonunknown);
		if ($taxonunknown eq "Unknown" || $taxonunknown eq "Environmental Unknown")
			{	print UNKNOWNS "$i\t";
				if (exists $ASV{$i}{'confidence'})
					{	print UNKNOWNS "$ASV{$i}{'confidence'}\t";
					}
				if (exists $ASV{$i}{'taxid'})
					{	print UNKNOWNS "$ASV{$i}{'taxid'}\t";
						my $string = $ASV{$i}{'taxid'};
						if ($string =~ m/\,/ || $string =~ m/\;/)
							{	my @unknown_tax_list;
								my @multi = split(',', $string);
								foreach my $j (@multi)
									{	$j =~ s/\ //;
										if ($j =~ m/\;/)
											{	my @multipli = split(';', $j);
												foreach my $k (@multipli)
													{	push(@unknown_tax_list, $k);
													}
											}
										else {push(@unknown_tax_list, $j);}
									}
								foreach my $entry (@unknown_tax_list)
									{	print UNKNOWNS "$TAXON{$entry}{'taxastring'}\t";
									}
								print UNKNOWNS "\n";
							}
						else
						{	print UNKNOWNS "$TAXON{$string}{'taxastring'}\n";
						}
					}
				else
					{	print UNKNOWNS "NO TAXID ASSIGNMENT\n";
					}
			}
	}
close(UNKNOWNS);

if ($options{d})
	{	system("mkdir IGNORING_ASVs");
		system("mv *IGNORE* IGNORING_ASVs/");
	}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
