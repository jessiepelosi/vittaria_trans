#!/usr/bin/perl -w

###########################################################################################################
# Program: runFisher_positive.pl
# Input: 1)MAPS output from actual data; 2)MAPS output from positive simulations;
# Output: summary of Fisher’s exact test
# Usage: perl runFisher_positive.pl
# Author: Barker Lab @ University of Arizona
# Last modified by Jessie Pelosi, Nov 9 2021
###########################################################################################################

#set the p-value with default 0.05
$bfp = 0.05;

%obs = ();
@mapsfiles = glob("*subtree.csv");
foreach $file (@mapsfiles)
{
        if ($file =~ m/(\S+)\_subtree\.csv/)
        {
                $tree = $1;
                $n = 0;
                print "$file \n";
                open FH1, '<', "$file";
            while (<FH1>)
            {
                    if (/^N\d+,\S+%,\S+%,(\d+),(\d+),\d+/)
                {
                                $b = $1;
                                $a = $2;
                                $n++;
                                $node = "N" . "$n";
                                $obs{$tree}{$node}{tl} = $a;
                                $obs{$tree}{$node}{bl} = $b;

                }
                }
                close FH1;
        }
}


$newFile = 0;
@csvfiles = glob("*.mapsMeanOut.csv");
foreach $file (@csvfiles)
{
    if ($file =~ m/(\S+)\.mapsMeanOut\.csv/)
    {
                $tree = $1;
                #$tree = uc($tree);
                $n = 0;
                print "$file \n";
                if (exists $obs{$tree}{"N1"}{tl})
                {
                        print "$file \n";
                        open FH1, '<', "$file";
                        while (<FH1>)
                        {
                        if (/^N\d+\,\S+\%\,\S+\%\,(\d+)\.\d+\,(\d+)\.\d+\,\S+\,\d+/)
                            {
                                        print "$file \n";
                                        $d = $1;
                                        $c = $2;
                                        $n++;
                                        $node = "N" . "$n";
                                        open OUT1,'>', "$tree.$n.fisher.R";
                                        print OUT1 "obs <- c($obs{$tree}{$node}{tl},$obs{$tree}{$node}{bl})\n";
                                        print OUT1 "null <- c($c,$d)\n";
                                        print OUT1 "mat <- cbind(obs,null)\n";
                                        print OUT1 "res <- fisher.test(mat,alternative=\"less\")\n";
                                        if ($newFile == 0)
                                        {
                                        print OUT1 "if (res\$p.value < $bfp){\nwrite(paste(\"$tree\\t$node\\t\",res\$p.value,\"\\t\",res\$estimate,\"\\t*\"),file=\"fisher.results.txt\")\n}\n";
                                            print OUT1 "if (res\$p.value >= $bfp){\nwrite(paste(\"$tree\\t$node\\t\",res\$p.value,\"\\t\",res\$estimate,\"\\tNS\"),file=\"fisher.results.txt\")\n}\n";
                                            $newFile = 1;
                                        }
                                        elsif ($newFile == 1)
                                        {
                                print OUT1 "if (res\$p.value < $bfp){\nwrite(paste(\"$tree\\t$node\\t\",res\$p.value,\"\\t\",res\$estimate,\"\\t*\"),file=\"fisher.results.txt\",append=TRUE)\n}\n";
                            print OUT1 "if (res\$p.value >= $bfp){\nwrite(paste(\"$tree\\t$node\\t\",res\$p.value,\"\\t\",res\$estimate,\"\\tNS\"),file=\"fisher.results.txt\",append=TRUE)\n}\n";
                                        }
                                        close OUT1;
                                        system "R CMD BATCH $tree.$n.fisher.R";
                                        unlink("$tree.$n.fisher.R");
                                        unlink("$tree.$n.fisher.Rout");
                            }
                        }
                }
    }
}
exit;
