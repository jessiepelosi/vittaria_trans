#!/usr/bin/perl -w

###########################################################################################################
# Program: simulateGeneTrees.3.0.pl
# Input: 1)An ultrametric species tree in newick format; 2) A sim.ctl file;
# Output: simulated gene trees
# Usage: perl simulateGeneTrees.3.0.pl
# Author: Barker Lab @ University of Arizona
# Last revised by Jessie Pelosi, University of Florida October 27, 2021
###########################################################################################################

# Parsing control file
open FH1, '<', "sim.ctl" or die("Missing sim.ctl: The control file must be named \"sim.ctl\"\n");
while (<FH1>)
{
        if (/^\w+/)
        {
                if (/treeFile\s+(\S+)/)
                {
                        $sTree = $1;
                }
                if (/nTrees\s+(\S+)/)
                {
                        $nTrees = $1;
                }
                if (/force_minimum\s+(\d+)/)
                {
                        $force_minimum = $1;
                }
                if (/fixed_lambda\s+(\S+)/)
                {
                        $fixed_lambda = $1;
                }
                if (/fixed_mu\s+(\S+)/)
                {
                        $fixed_mu = $1;
                }
                if (/outFile\s+(\S+)/)
                {
                        $outFile= $1;
                }
                if (/bootStrap\s+(\d+)/)
                {
                        $bootStrap= $1;
                }
                if (/nReps\s+(\d+)/)
                {
                        $nReps= $1;
                }
                if (/root_distribution\s+(\S+)/)
                {
                        $root_distribution = $1;
                        $geom_p = 1/$root_distribution;
                }
                if (/wgd_retention_rate\s+(.+)/)
                {
                        $rateString = $1;
                        @wgd_retention_rate = ();
                        @wgd_retention_rate = split(/ / ,$rateString);
                        for $i (0..(scalar(@wgd_retention_rate) - 1))
                        {
                                print "WGD $i has retention rate $wgd_retention_rate[$i]\n";
                        }
                }
                if (/wgd_time_before_divergence\s+(.+)/)
                {
                        $timeString = $1;
                        @wgd_time_before_divergence = ();
                        @wgd_time_before_divergence = split(/ / ,$timeString);
                }
        }
}
close FH1;

##########################################################################################
#Get list of all species in tree
%sList = ();
$nTax = 0;
open FH1, '<', "$sTree";
while (<FH1>)
{
        $treeString = $_;
        chomp $treeString;
        $treeString_copy = $treeString;
        $treeString_wgd = $treeString;
        $treeString_copy =~ s/\;//;
        $treeString_copy =~ s/WGD\d+//g;
        $treeString =~ s/WGD\d+//g;
        $treeString =~ s/\:\d+\.\d+/,/g;
        $treeString =~ s/\)/,/g;
        $treeString =~ s/\(/,/g;
        $treeString =~ s/\;/,/g;
        @temp = ();
        @temp = split(/\,/, $treeString);
        foreach $chunk (@temp)
        {
                if ($chunk =~ m/(\S+)/ && $chunk !~ m/\,/)
                {
                        $species = $1;
                        if (! exists $sList{$species})
                        {
                                $sList{$species} = 1;
                                $nTax++;
                        }
                }
        }
}
close FH1;
##########################################################################################

##########################################################################################
#Prepare a discrete distribution of gene count probability breakpoints for the root distribution up to floating point precision
if ($root_distribution > 1)
{
        $still_floating = 1;
        $k = 1;
        @prob_rootCount = ();
        $cum_prob = 0;
        while ($still_floating == 1)
        {
                $prob = ((1- $geom_p) ** ($k - 1)) * $geom_p;
                $cum_prob = $cum_prob + $prob;
                if ($prob >= 1e-8)
                {
                        push @prob_rootCount, $cum_prob;
                }
                elsif ($prob < 1e-8)
                {
                        push @prob_rootCount, 1;
                        $still_floating = 0;
                }
                $k++;
        }
        print "Cumulative probability distribution of gene copies at root\n@prob_rootCount\n";
}
elsif($root_distribution == 1)
{
        print "Assuming a single gene copy across all gene families at the root of the species tree\n";
}
elsif($root_distribution < 1)
{
        die ("Distribution of gene family sizes at root of species tree must have a geometric mean of 1 or greater ... exiting\n");
}
##########################################################################################

##########################################################################################
#Find the node to undergo a WGD
#Include the branch length before the node too

@wgd_subtree_array = ();
@wgd_taxa_array = ();
for $w (0..(scalar(@wgd_retention_rate) - 1))
{
        $wgd_subtree = "";
        $wgd_node_bl = "";
        $bl_switch = 0;
        %wgd_taxa = ();
        $check = 0;
        @temp = ();
        @temp = split(// ,$treeString_wgd);
        for $i (0..(scalar(@temp)-1))
        {
                if ($temp[$i] eq "W" and $temp[($i+1)] eq "G" and $temp[($i+2)] eq "D" and $temp[($i+3)] eq "$w")
                {
                        $pc = 0;
                        $po = 0;
                        for $i2 (0..($i - 1))
                        {
                                if ($temp[$i  - 1 - $i2] eq ")")
                                {
                                        $po++;
                                }
                                elsif ($temp[$i - 1 - $i2] eq "(")
                                {
                                        $pc++;
                                }
                                if ($check == 0)
                                {
                                        $wgd_subtree = $temp[$i - 1 - $i2] . $wgd_subtree;
                                }
                                if ($pc == $po)
                                {
                                        $check = 1;
                                }
                        }
                }
                if ($bl_switch == 1)
                {
                        if ($temp[$i] eq ")" or $temp[$i] eq "(" or $temp[$i] eq ",")
                        {
                                $bl_switch = 0;
                        }
                        else
                        {
                                $wgd_node_bl = $wgd_node_bl . $temp[$i];
                        }
                }
                if ($temp[$i] eq "$w" and $temp[($i-1)] eq "D" and $temp[($i-2)] eq "G" and $temp[($i-3)] eq "W")
                {
                        $bl_switch = 1;
                }
        }
        $wgd_subtree = $wgd_subtree . $wgd_node_bl;
        $tempTree = $wgd_subtree;
        $tempTree =~ s/WGD\d+//g;
        $tempTree =~ s/:\d+\.\d+//g;
        $tempTree =~ s/[\(\)]//g;
        @temp = ();
        @temp = split(/,/ ,$tempTree);
        foreach $temp_tax (@temp)
        {
                $wgd_taxa{$temp_tax} = 1;
                $wgd_taxa_array[$w]{$temp_tax} = 1;
        }
        $wgd_subtree =~ s/WGD\d+//g;
        print "WGD $w will duplicate subtree: $wgd_subtree at $wgd_time_before_divergence[$w] before MRCA\n";
        $wgd_subtree_array[$w][0] = $wgd_subtree;
        $wgd_subtree_array[$w][1] = $wgd_subtree;
}


##########################################################################################

##########################################################################################
#All topologies will go to outFile to be analyzed directly with MAPS
open OUT1, '>', "$outFile";

$wgd_switch = 0;

$i = 1;
$wgd_i = 1;
while ($i <= $nTrees)
{
    if ($root_distribution > 1)
    {
        $sampletree = getStreeRoot();
        }
    elsif ($root_distribution == 1)
        {
                $sampletree = "$treeString_copy;";
        }
        #####
        #Check that retention rates are non-zero. If so sample a species tree with WGD
        for $z (0..(scalar(@wgd_retention_rate) - 1))
        {
                if ($wgd_retention_rate[$z] > 0)
                {
                        $wgd_switch = 1;
                }
        }
        if ($wgd_switch == 1)
        {
                $sampletree = getWGDTree($sampletree);
                print "$sampletree \n";
        }
        #####
        if ($force_minimum == 0)
        {
                system "java -jar jprime-0.3.6.jar GuestTreeGen -n 1 -max 1000 -maxper 100 \"$sampletree\" $fixed_lambda $fixed_mu 0.0 fam.$i";
        }
        elsif ($force_minimum == 1)
        {
                system "java -jar jprime-0.3.6.jar GuestTreeGen -n 1 -max 1000 -minper 1 -maxper 100 \"$sampletree\" $fixed_lambda $fixed_mu 0.0 fam.$i";
        }
        system "java -jar jprime-0.3.6.jar BranchRelaxer -o fam.$i.relaxed.tre fam.$i.pruned.tree CONSTANT 1";
        open FH1, '<', "fam.$i.pruned.leafmap";
        open FH2, '<', "fam.$i.relaxed.tre";
        while (<FH2>)
        {
                $treestring = $_;
            chomp $treestring;
        }
        close FH2;
        while (<FH1>)
        {
                if (/(G(\d+)\_\d+)\s+(\S+)/)
            {
                $leaf = $1;
                $tipnum = $2;
                $tax = $3;
                $treestring =~ s/$leaf/$tax\_ID_$tipnum/;
            }
        }
        close FH1;
#Strip branch lengths from tree. Remove unique identifiers for MAPS.
        $treestring =~ s/\:\d+\.\d+//g;
        $treestring =~ s/E-8//g;
        $treestring =~ s/\_ID\_\d+//g;
        $treestring =~ s/\_\d+//g;
#Check tree has at least on leaf per species
        $treeCopy = $treestring;
        $treeCopy =~ s/\)/,/g;
        $treeCopy =~ s/\(/,/g;
        $treeCopy =~ s/\;/,/g;
        $theseTax = 0;
        %temp_sList = ();
        @temp = ();
        @temp = split(/\,/, $treeCopy);
        foreach $chunk (@temp)
        {
                if ($chunk =~ m/(\S+)/ && $chunk !~ m/\,/)
                {
                        $species = $1;
                        if (! exists $temp_sList{$species} && exists $sList{$species})
                        {
                                $temp_sList{$species} = 1;
                                $theseTax++;
                        }
                }
        }
#Clean up temporary files
        if ($theseTax == $nTax)
        {
                print OUT1 "$treestring\n";
                unlink glob("fam.$i.unpruned.*");
                unlink glob("fam.$i.pruned.*");
                unlink glob("fam.$i.relaxed.tre*");
                $i++;
                $treestring = "";
        }
        elsif ($theseTax != $nTax)
        {
                unlink glob("fam.$i.unpruned.*");
                unlink glob("fam.$i.pruned.*");
                unlink glob("fam.$i.relaxed.tre*");
                $treestring = "";
        }
}

close OUT1;
##########################################################################################

##########################################################################################
#Resample with replacement from set of simulated topologies to generate bootstraps
#The names of bootstrapped dataset files are fixed
if ($bootStrap == 1)
{
        @simulatedTrees = ();
        open FH1, '<', "$outFile";
        {
                while (<FH1>)
                {
                        $t = $_;
                        chomp $t;
                        push @simulatedTrees, $t;
                }
        }
        close FH1;
        for $i (1..$nReps)
        {
                open OUT1, '>', "boot.$i.txt";
                for (1..$nTrees)
                {
                        $r = int(rand($nTrees - 1));
                        $sample = $simulatedTrees[$r];
                        print OUT1 "$sample\n";
                }
                close OUT1;
        }
}
##########################################################################################

##########################################################################################
#Sample a gene tree where the number of lineages at the root are drawn from a geometric distribution with a given mean
sub getStreeRoot
{
    $root_prob = rand(1);
    $found_k = 0;
        $this_k = 1;
    OUTER: while ($found_k == 0)
        {
                for $k (1..(scalar(@prob_rootCount) -1))
        {
                if ($root_prob <= $prob_rootCount[$k])
                {
                        $this_k = $k;
                        $found_k = 1;
                        last OUTER;
                }
        }
        }
        if ($this_k > 1)
        {
                $tree_i = "";
                for $l (1..($this_k))
                {
                        $thisTree = $treeString_copy;
                        foreach $leaf (keys %sList)
                        {
                                $thisTree =~ s/$leaf/$leaf\_$l/;
                        }
                        if ($l == 1)
                        {
                                $tree_i = "$thisTree" . ":0.00000001";
                        }
                        elsif ($l == 2)
                        {
                                $tree_i = "$tree_i" . "," . "$thisTree" . ":0.00000001";
                        }
                        elsif ($l > 2)
                        {
                                $tree_i = "(" . "$tree_i" . ")" . ":0.00000001" . "," . "$thisTree" . ":0.00000001";
                        }
                }
                $tree_i = "(" . "$tree_i" . ")" . ";";
        }
        elsif ($this_k == 1)
        {
                $tree_i = "$treeString_copy" . ";";
        }
        return $tree_i;
}
##########################################################################################
#Embed WGDs into the species tree
#IF there are multiple copies at the root, drop those appended ids. Individual copies will get new ones.
sub getWGDTree
{
        $wgdTree = $_[0];
        $wgdTree =~ s/\_\d+//g;
        for $w (0..(scalar(@wgd_retention_rate) - 1))
        {
                $p_wgd = rand(1);
                if ($p_wgd <= $wgd_retention_rate[$w])
                {
                        $wgd_duplicated_subtree = "";
                        $left = $wgd_subtree_array[$w][0];
                        $right = $wgd_subtree_array[$w][1];
                        if ($left =~ m/\S+\)\:(\d+\.\d+)$/)
                        {
                                $old_bl = $1;
                                $new_bl = $old_bl - $wgd_time_before_divergence[$w];
                        }
                        $left =~ s/$old_bl/$wgd_time_before_divergence[$w]/;
                        $right =~ s/$old_bl/$wgd_time_before_divergence[$w]/;
                        $wgd_duplicated_subtree = "($left,$right):$new_bl";

                        $temp_wgd_subtree = $wgd_subtree_array[$w][0];
                        $temp_wgd_subtree =~ s/\(/\\(/g;
                        $temp_wgd_subtree =~ s/\)/\\)/g;
                        $temp_wgd_duplicated_subtree = $wgd_duplicated_subtree;
                        $temp_wgd_duplicated_subtree =~ s/\(/\\(/g;
                        $temp_wgd_duplicated_subtree =~ s/\)/\\)/g;
                        $wgdTree =~ s/$temp_wgd_subtree/$temp_wgd_duplicated_subtree/g;
                        $wgdTree =~ s/\\//g;
                        $wgdTree = $wgdTree;
                }
        }
        #Add unique ids to end of each leaf so jprime is cool with it
        $id = 1;
        $maybeLeaf = "";
        $wgdTree_new = "";
        @finalTreeArray = split(//,$wgdTree);
        for $this_i (0..(scalar(@finalTreeArray)-1))
        {
                if ($finalTreeArray[$this_i] !~ m/[\(\)\,\;\:]/)
                {
                        $maybeLeaf = $maybeLeaf . "$finalTreeArray[$this_i]";
                }
                elsif ($finalTreeArray[$this_i] =~ m/[\(\)\,\;\:]/)
                {
                        if ($maybeLeaf !~ m/\d+\.\d+/)
                        {
                                if ($maybeLeaf =~ m/\S+/)
                                {
                                        $wgdTree_new = $wgdTree_new . "_" . "$id";
                                        $id++;
                                        $maybeLeaf = "";
                                }
                        }
                        elsif ($maybeLeaf =~ m/\d+\.\d+/)
                        {
                                $maybeLeaf = "";
                        }
                }
                $wgdTree_new = $wgdTree_new . "$finalTreeArray[$this_i]";
        }
        return $wgdTree_new;
}
##########################################################################################
exit;
