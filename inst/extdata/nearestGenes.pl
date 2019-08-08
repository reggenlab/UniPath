#!/usr/bin/perl -w

$#ARGV == 2 || die "Two arguments required: filename_of_coordinates and refseq_file_name and output_filename ";

$larg_dist = 1e10;

# print $ARGV[0] . " 0 " . $ARGV[1] . " 1 " . $ARGV[2] . "\n" ;
 
$genfile = $ARGV[1] ; 
 open( genf1, "$genfile") || die "couldn't open refeq file\n";
@gen1 = <genf1>;
close(genf1);

open( outf1 , ">tempfile1") ;
for($i = 0 ;  $i < scalar(@gen1) ; $i++)
{
$temps = $gen1[$i] ;
chomp($temps) ;
@ss = split('\t' , $temps) ;
print outf1 "$ss[2]\t$ss[4]\t$ss[5]\t$ss[3]\t$ss[1]\t$ss[12]\t$ss[8]\n" ; 
}
close(outf1) ;




open(GENE, "tempfile1") || die "couldn't open refeq file\n";;
@gen = <GENE>;
close(GENE);

open( REG , $ARGV[0]) ;
@reg = <REG> ; 
close(REG) ;

%regs_by_chrom = ();
%nreg_by_chrom = ();

$r = 0;
#print $reg[1]. "\n" .$reg[2] ;

foreach (@reg)
{
   if ( /^chr/ )
    {
	@seg = split;
	++$nreg_by_chrom{$chrom = $seg[0]};
	$regs_by_chrom{$chrom}[$nreg_by_chrom{$chrom} - 1] = [@seg, $r++];
    }
}

%genes_by_chrom = ();
%ngene_by_chrom = ();

$g = 0;
foreach (@gen)
{
    if ( /^chr/ )
    {
	@seg = split;
	++$ngene_by_chrom{$chrom = $seg[0]};
	$genes_by_chrom{$chrom}[$ngene_by_chrom{$chrom} - 1] = [@seg, $g++];
    }
}

@closest[0..(@reg-1)] = (0);
foreach $chrom (keys %regs_by_chrom)
{

    @regs = @{$regs_by_chrom{$chrom}};
    $nreg = $nreg_by_chrom{$chrom};

    @genes = @{$genes_by_chrom{$chrom}};
    $ngene = $ngene_by_chrom{$chrom};

    print STDERR "Starting on $chrom\n";


    for ($i = 0; $i < $nreg; $i++)
    {
	$cdL = -$larg_dist;   # left
	$cdR = $larg_dist;    # right
	$cdI = $larg_dist;    # intronic
	($cgL, $cgR, $cgI) = ();

	for ($j = 0; $j < $ngene; $j++)
	{
	    if ( is_intronic($regs[$i], $genes[$j]) )
	    {
		$d = abs(reg_gen_dist($regs[$i], $genes[$j]));
		if ($d < $cdI)
		{
		    $cdI = $d;
		    $cgI = $genes[$j][5];
		}
	    }
	    else
	    {
		$d = reg_gen_dist($regs[$i], $genes[$j]);
		if ($d <= 0)
		{		    
		    if ($d > $cdL)
		    {
			$cdL = $d;
			$cgL = $genes[$j][5];
		    }
		}
		else
		{		    
		    if ($d < $cdR)
		    {
			$cdR = $d;
			$cgR = $genes[$j][5];
		    }
		}
	    }
	}

	@cdists = (abs($cdL), $cdR, $cdI);
	@cgenes = ($cgL, $cgR, $cgI);

# Choose the two closest genes from among the left-gene, the right-gene and the intron-gene (if it exists).
	@cgen_idx = (sort {$cdists[$a] <=> $cdists[$b]} (0,1,2));
	@cg       = @cgenes[@cgen_idx[0,1]];
        @cd       = @cdists[@cgen_idx[0,1]];

	if (defined $cg[1])
	{
# if there are multiple isoforms, it could happen that, say, $cgR = $cgI
#	    $cg[1] eq $cg[0] && defined ($tmp_g = $cgenes[$cgen_idx[2]]) &&
#		($cg[1] = $tmp_g);
            $cg[1] eq $cg[0] && defined ($tmp_g = $cgenes[$cgen_idx[2]]) &&
		($cg[1] = $tmp_g) && ($cd[1] = $cdists[$cgen_idx[2]]);
	}
	else
	{
# a telomeric CNS would have only a right neighbour, or only a left neighbour
#	    $cg[1] = $cg[0];
            $cg[1] = $cg[0];
            $cd[1] = $cd[0];
	}

#	$closest[$regs[$i][-1]] = [@cg]; #last element of reg. is the index
        $closest[$regs[$i][-1]] = [$cg[0],$cd[0],$cg[1],$cd[1]]; #last element of reg. is the index
    }

}

$i = 1;
open(outf1 , ">$ARGV[2]") ;
foreach $c_ref (@closest)
{
    print outf1 join("\t", @$c_ref), "\n";
$i = $i + 1 ;
}


sub reg_gen_dist
{

    my ($reg_ref, $gen_ref) = @_;
    my ($rs, $re, $gs, $ge, $dmin);

    $rs = $reg_ref->[1] + 1;
    $re = $reg_ref->[2];
    $gs = $gen_ref->[1] + 1;
    $ge = $gen_ref->[2];

#    return (sort {abs($a) <=> abs($b)} ($gs-$rs, $gs-$re, $ge-$rs, $ge-$re))[0];
    
    ($gen_ref->[3] eq '+') && return (sort {$a <=> $b} ($gs-$rs, $gs-$re))[0];
    ($gen_ref->[3] eq '-') && return (sort {$a <=> $b} ($ge-$rs, $ge-$re))[0];
}


sub is_intronic
{

    my ($reg_ref, $gen_ref) = @_;

    my ($rs, $re, $gs, $ge);

    $rs = $reg_ref->[1] + 1;
    $re = $reg_ref->[2];
    $gs = $gen_ref->[1] + 1;
    $ge = $gen_ref->[2];

    $rs > $ge || $gs > $re || return 1;
    return 0;

}
