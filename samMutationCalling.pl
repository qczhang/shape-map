#! /usr/bin/perl
#
use strict;
use warnings;

my $_debug = 1;

use Data::Dumper;

my $inputSam = shift;
my $outputMut = shift;
my $outputStat = shift;
my $fasta = shift;
my $icSHAPE = shift;

my %seq_fasta = (); 
my %seq_icSHAPE = ();
my %seq_mut = ();
my %seqMut_stat = (); 
my %allMut_stat = ();

&main ( $inputSam, $outputMut, $outputStat, icSHAPE => $icSHAPE, fasta => $fasta );

sub main
{
    &init ();

    my $samFileList = shift;
    my $mutFile = shift;
    my $statFile = shift;
    my %parameters = @_;

    &readFasta ( $parameters{fasta}, simple => 1 );
    &readIcSHAPE ( $parameters{icSHAPE}, ) if ( defined $parameters{icSHAPE} );

    my @samFiles = split ( /:/, $samFileList );
    foreach my $samFile ( @samFiles ) { &readSam ( $samFile ); }

    &statAndPrint ( $mutFile, icSHAPE => $icSHAPE, lowCut => 0.05, highCut => 0.4 );
    &outputMutStat ( $statFile, icSHAPE => $icSHAPE );
}

sub init
{
    if ( not $inputSam ) { die "Usage: $0 input_sam output_stat fasta shape\n"; }
    if ( ( not defined $outputMut ) or ( $outputMut eq "NULL" ) ) {  $outputMut = "output.stat";  }
    if ( ( not defined $fasta ) or ( $fasta eq "NULL" ) ) { $fasta = "/home/qczhang/database/ensembl/current/mouse/gtf/transcriptome.fa"; }
    if ( ( not defined $icSHAPE ) or ( $icSHAPE eq "NULL" ) ) { $icSHAPE = "/home/qczhang/shape-seq/new/analysis/all.polya/LIB_NAI-LIB_DMSO.PolyA.invitro.valid.enrich"; }
}

sub readSam
{
    my $samFile = shift;

    my $lineCount = 0;
    open ( SAM, $samFile ) or die ( "Error in reading sam file $samFile!\n" );
    print STDERR "read sam file $samFile...\n";
    while ( my $line = <SAM> ) {
        next if ( $line =~ /^#/ );
        if ( $line =~ /^@/ ) { }
        else {
            $lineCount++; print STDERR "line: $lineCount\n\t", `date` if ( $lineCount % 1000000 == 0 );

            my @data = split ( /\t/, $line );
            my $tag = $data[1];
            if ( ( not $tag ) or ( $tag == 99 ) ) { ## so far we only use read1
                my $mdString = ""; 
                for ( my $idx = 11; $idx < scalar ( @data ); $idx++ ) { if ( $data[$idx] =~ /MD:/ ) { $mdString = substr ( $data[$idx], 5 ); last; } }
                &parseMut ( $data[2], $data[3], $data[5], $data[9], $mdString );
            }
            elsif ( ( $tag == 99 ) or ( $tag == 147 ) )  {
                next;

                $line = <SAM>;
                my @data2 = split ( /\t/, $line );
                if ( $data[0] ne $data2[0] ) { print STDERR "ERROR! consecutive reads not the same...skipped. Sort SAM files first?\n"; next; }
                if ( $data[2] ne $data2[2] ) { print STDERR "ERROR! consecutive reads not mapped to the same targets...skipped.\n"; next; }

                my $mdString = ""; my $mdString2 = "";
                for ( my $idx = 11; $idx < scalar ( @data ); $idx++ ) { if ( $data[$idx] =~ /MD:/ ) { $mdString = substr ( $data[$idx], 5 ); last; } } 
                for ( my $idx = 11; $idx < scalar ( @data2 ); $idx++ ) { if ( $data2[$idx] =~ /MD:/ ) { $mdString2 = substr ( $data2[$idx], 5 ); last; } }
                &parseMut2 ( $data[2], $data[3], $data2[3], $data[5], $data2[5], $data[9], $data2[9], $mdString, $mdString2 );
            }
        }
    }
    close SAM;

    return $lineCount;
}

sub parseMut
{
    my $seqID = shift; my $pos = shift; my $cigar = shift; my $readSeq = shift; my $md = shift;

    my ( $ref_match, $ref_matchSize ) = _parseCigar ( $cigar );
    my ( $ref_op, $ref_opSize ) = _parseMD ( $md );

    my ( $headSoftClip, $ref_insertions, $ref_insertionSize ) = collectAlignInfoCIGAR ( $seqID, $pos, $cigar, $readSeq, $ref_match, $ref_matchSize );
    collectAlignInfoMD ( $seqID, $pos, $readSeq, $headSoftClip, $ref_op, $ref_opSize, $ref_insertions, $ref_insertionSize );

    1;
}

sub collectAlignInfoCIGAR
{
    my $seqID = shift; my $pos = shift; my $cigar = shift; my $readSeq = shift;
    my $ref_match = shift; my $ref_matchSize = shift;

    my $refPos = $pos - 1;      ## read alignment position ( pos ) is 1-indexed, but refPos is 0-indexed
    my $readPos = 0;            ## read position is 0-indexed on read ( column 10 )
    my $headSoftClip = 0;       ## whether the leading is softclipped
    my @insertions = (); my @insertionSize = ();
    for ( my $idx = 0; $idx < scalar ( @{$ref_match} ); $idx++ ) {
        if ( $ref_match->[$idx] eq "S" )  { 
            if ( $idx == 0 ) { 
                $seq_mut{$seqID}{"sclip"}[$refPos] = _append ( $seq_mut{$seqID}{"sclip"}[$refPos], substr ( $readSeq, $readPos, $ref_matchSize->[$idx] ) . "-" );
                $headSoftClip = $ref_matchSize->[$idx];
            }
            elsif ( $idx == scalar( @{$ref_match} -1 ) ) { 
                $seq_mut{$seqID}{"sclip"}[$refPos] = _append ( $seq_mut{$seqID}{"sclip"}[$refPos], "-" . substr ( $readSeq, $readPos, $ref_matchSize->[$idx] ) );
            }
            else { print STDERR "Warning! clipping inside an alignment!\n"; print STDERR "\t$cigar\t$readSeq\n"; }
            $readPos += $ref_matchSize->[$idx];
        }
        if ( $ref_match->[$idx] eq "H" )  { 
            if ( $idx == 0 ) { $seq_mut{$seqID}{"hclip"}[$refPos] = _append ( $seq_mut{$seqID}{"hclip"}[$refPos], "csna-" ); }
            elsif ( $idx == scalar( @{$ref_match} -1 ) ) { $seq_mut{$seqID}{"hclip"}[$refPos] = _append ( $seq_mut{$seqID}{"hclip"}[$refPos], "-csna" ); }
            else { print STDERR "Warning! clipping inside an alignment!\n"; print STDERR "\t$cigar\t$readSeq\n"; }
        }
        elsif ( $ref_match->[$idx] eq "I" ) {
            $seq_mut{$seqID}{"insertion"}[$refPos] = _append ( $seq_mut{$seqID}{"insertion"}[$refPos], substr ( $readSeq, $readPos, $ref_matchSize->[$idx] ) );
            $readPos += $ref_matchSize->[$idx];
            push ( @insertions, $refPos );
            push ( @insertionSize, $ref_matchSize->[$idx] );
        }
        elsif ( $ref_match->[$idx] eq "P" ) { }
        elsif ( $ref_match->[$idx] eq "D" ) {
            $refPos += $ref_matchSize->[$idx];     
            # $seq_mut{$seqID}{"deletion"}[$refPos] = _append ( $seq_mut{$seqID}{"deletion"}[$refPos], $ref_matchSize->[$idx] ); # subject to realignment
        }
        elsif ( $ref_match->[$idx] eq "X" ) {
            $seq_mut{$seqID}{"mismatch"}[$refPos] = _append ( $seq_mut{$seqID}{"mismatch"}[$refPos], substr ( $readSeq, $readPos, $ref_matchSize->[$idx] ) ); 
            ## not tested
            $refPos += $ref_matchSize->[$idx];     
            $readPos += $ref_matchSize->[$idx];
        }
        elsif ( ( $ref_match->[$idx] eq "=" ) or ( $ref_match->[$idx] eq "M" ) ) {
            $refPos += $ref_matchSize->[$idx];     
            $readPos += $ref_matchSize->[$idx];
        }
    }

    return ( $headSoftClip, \@insertions, \@insertionSize );
}

sub collectAlignInfoMD
{
    my $seqID = shift; my $pos = shift; my $readSeq = shift; my $headSoftClip = shift;
    my $ref_op = shift; my $ref_opSize = shift;
    my $ref_insertions = shift; my $ref_insertionSize = shift;

    my $refPos = $pos - 1;               ## point to 0-base index 
    my $readPos = $headSoftClip;         ## the same

    if ( $ref_op->[0] ) {
        if ( $ref_op->[0] =~ /^\^/ ) {
            my $delLen = length ( $ref_op->[0] ) - 1;
            for ( my $idxPos = 0; $idxPos < $delLen-1; $idxPos++ ) { $seq_mut{$seqID}{"same"}[$refPos+$idxPos]++; }
            $refPos += $delLen;
            if ( ( defined $ref_op->[1] ) and ( not defined $ref_opSize->[1] ) ) { print STDERR "error!\n"; return -1; }
            if ( ( defined $ref_op->[1] ) and ( $ref_op->[1] ) and ( not $ref_opSize->[0] ) ) { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            else {
                my $matchInLocal = _localReAlignment ( $seqID, $refPos-$delLen, $ref_op->[0] );
                print STDERR "multiple matching? - $matchInLocal\n" if ( $_debug );
                if ( $matchInLocal == 1 ) {
                    print STDERR "\tmatch only once, no realignment needed.\n" if ( $_debug );
                    $seq_mut{$seqID}{"deletion"}[$refPos-1] = _append ( $seq_mut{$seqID}{"deletion"}[$refPos-1], $ref_op->[0] );
                    ## you want to label the deletion even at the last base of deleted fragment
                }
                else { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            }
        }
        elsif ( $ref_op->[0] =~ /[A-Z]/ ) {
            $refPos++; $readPos++;
            if ( ( defined $ref_op->[1] ) and ( not defined $ref_opSize->[1] ) ) { print STDERR "error!\n"; return -1; }
            if ( ( defined $ref_op->[1] ) and ( $ref_op->[1] ) and ( not $ref_opSize->[1] ) ) { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            else { $seq_mut{$seqID}{"mutation"}[$refPos-1] = _append ( $seq_mut{$seqID}{"mutation"}[$refPos-1], substr ( $readSeq, $readPos-1, 1 ) ); }
        }

        $ref_op->[0] = "";
        shift @{$ref_opSize};
    }

    for ( my $idx = 0; $idx <= scalar ( @{$ref_op} ); $idx++ ) {
        if ( ( not defined $ref_op->[$idx] ) or ( not $ref_op->[$idx] ) ) {
        }
        elsif ( $ref_op->[$idx] =~ /^\^/ ) {
            my $delLen = length ( $ref_op->[$idx] ) - 1;
            for ( my $idxPos = 0; $idxPos < $delLen-1; $idxPos++ ) { $seq_mut{$seqID}{"same"}[$refPos+$idxPos]++; } #deletion only record at the last position
            $refPos += $delLen;
            if ( ( defined $ref_op->[$idx+1] ) and ( $ref_op->[$idx+1] ) and ( not $ref_opSize->[$idx] ) ) { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            else {
                my $matchInLocal = _localReAlignment ( $seqID, $refPos-$delLen, $ref_op->[$idx] );
                print STDERR "multiple matching? - $matchInLocal\n" if ( $_debug );
                if ( $matchInLocal == 1 ) {
                    print STDERR "\tmatch only once, no realignment needed.\n" if ( $_debug );
                    $seq_mut{$seqID}{"deletion"}[$refPos-1] = _append ( $seq_mut{$seqID}{"deletion"}[$refPos-1], $ref_op->[$idx] );
                    ## you want to label the deletion even at the last base of deleted fragment
                }
                else { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            }
        }
        elsif ( $ref_op->[$idx] =~ /[A-Z]/ ) {
            $refPos++; $readPos++;
            if ( ( defined $ref_op->[$idx+1] ) and ( $ref_op->[$idx+1] ) and ( not $ref_opSize->[$idx] ) ) { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            else { $seq_mut{$seqID}{"mutation"}[$refPos-1] = _append ( $seq_mut{$seqID}{"mutation"}[$refPos-1], substr ( $readSeq, $readPos-1, 1 ) ); }
        }

        if ( $ref_opSize->[$idx] ) {
            for ( my $idxPos = 0; $idxPos < $ref_opSize->[$idx]; $idxPos++ ) { $seq_mut{$seqID}{"same"}[$refPos+$idxPos]++; }
            $refPos += $ref_opSize->[$idx];
            $readPos += $ref_opSize->[$idx];

            while ( ( scalar (@{$ref_insertions} ) ) and ( $refPos >= $ref_insertions->[0] ) )  {
                $readPos += $ref_insertionSize->[0];
                shift ( @{$ref_insertions} );
                shift ( @{$ref_insertionSize} );
            }
        }
    }

    1;
}
sub _parseMut2
{
}

sub statAndPrint
{
    my $outFile = shift;
    my %parameters = @_;

    open ( OUT, ">$outFile" );
    print OUT "seqID\tindex\tbase\ticSHAPE\tmutFqA\tmutFqT\tmutFqG\tmutFqC\tmutFqDel\tmutFqTotal\tnoMut\ttotal\tmutPbA\tmutPbT\tmutPbG\tmutPbC\tmutPbDel\tmutPbTotal\n";
    foreach my $seqID ( sort {$a cmp $b} ( keys %seq_mut ) ) {
        if ( not defined $seq_fasta{$seqID} ) { print STDERR "ERROR! Sequence of $seqID not found. ...skipped\n"; next; }

        print OUT "## -- $seqID --\n";
        my $fasta = uc ( $seq_fasta{$seqID} );
        for ( my $idx = 0; $idx < length ( $fasta ); $idx++ ) {
            my $base = substr ( $fasta, $idx, 1 );
            my $icSHAPE = ( defined $seq_icSHAPE{$seqID} ) ? $seq_icSHAPE{$seqID}[$idx] : "NULL"; 
            print OUT $seqID, "\t", $idx+1, "\t", $base, "\t", $icSHAPE;

            my $total_mut = 0; my $total = 0;
            my $insA = 0; my $insT = 0; my $insG = 0; my $insC = 0;
            if ( defined $seq_mut{$seqID}{mutation}[$idx] ) {
                $insA = () = $seq_mut{$seqID}{mutation}[$idx] =~ /A/gi;
                $insT = () = $seq_mut{$seqID}{mutation}[$idx] =~ /T/gi;
                $insG = () = $seq_mut{$seqID}{mutation}[$idx] =~ /G/gi;
                $insC = () = $seq_mut{$seqID}{mutation}[$idx] =~ /C/gi;
                updateSeqMut ( $seqID, $base, "A", $insA, $icSHAPE, lowCut => $parameters{lowCut}, highCut => $parameters{highCut} );
                updateSeqMut ( $seqID, $base, "T", $insT, $icSHAPE, lowCut => $parameters{lowCut}, highCut => $parameters{highCut} );
                updateSeqMut ( $seqID, $base, "G", $insG, $icSHAPE, lowCut => $parameters{lowCut}, highCut => $parameters{highCut} );
                updateSeqMut ( $seqID, $base, "C", $insC, $icSHAPE, lowCut => $parameters{lowCut}, highCut => $parameters{highCut} );
            }

            my $insD = 0;
            if ( defined $seq_mut{$seqID}{deletion}[$idx] ) {
                $insD = () = $seq_mut{$seqID}{deletion}[$idx] =~ /;/gi; $insD++;
                updateSeqMut ( $seqID, $base, "deletion", $insD, $icSHAPE, lowCut => $parameters{lowCut}, highCut => $parameters{highCut} );
            }

            my $same = 0; 
            if ( defined $seq_mut{$seqID}{same}[$idx] ) {  
                $same = $seq_mut{$seqID}{same}[$idx];
                updateSeqMut ( $seqID, $base, $base, $same, $icSHAPE, lowCut => $parameters{lowCut}, highCut => $parameters{highCut} ); 
            }

            $total_mut += $insA + $insT + $insG + $insC + $insD;
            $total = $same + $total_mut;
            print OUT "\t", $insA, "\t", $insT, "\t", $insG, "\t", $insC, "\t", $insD, "\t", $total_mut, "\t", $same, "\t", $total;
            if ( $total ) 
                {  print OUT "\t", sprintf ( "%.4f", $insA/$total ), "\t", sprintf ( "%.4f", $insT/$total ), "\t", sprintf ( "%.4f", $insG/$total ), "\t", sprintf ( "%.4f", $insC/$total ), "\t", sprintf ( "%.4f", $insD/$total ), "\t", sprintf ( "%.4f", $total_mut/$total ), "\n"; }
            else { print OUT "\t-\t-\t-\t-\t-\t-\n"; }
        }
    }
    close OUT;

    1;
}

sub updateSeqMut
{
    my $seqID = shift; my $base = shift; my $mutType = shift; my $count = shift; my $icSHAPE = shift;
    my %parameters = @_;

    if ( not defined $seqMut_stat{$seqID}{$base}{$mutType} )  { $seqMut_stat{$seqID}{$base}{$mutType} = $count; }
    else  { $seqMut_stat{$seqID}{$base}{$mutType} += $count; }
    if ( not defined $allMut_stat{$base}{$mutType} )  { $allMut_stat{$base}{$mutType} = $count; }
    else  { $allMut_stat{$base}{$mutType} += $count; }

    if ( $icSHAPE ne "NULL" ) {
        my $cut = 0;
        if ( $icSHAPE > $parameters{highCut} ) { $cut = "highCut"; }
        elsif ( $icSHAPE < $parameters{lowCut} ) { $cut = "lowCut"; }

        if ( $cut ) {
            if ( not defined $seqMut_stat{$seqID}{$base}{$cut}{$mutType} )  { $seqMut_stat{$seqID}{$base}{$cut}{$base} = $count; }
            else  { $seqMut_stat{$seqID}{$base}{$cut}{$mutType} += $count; }
            if ( not defined $allMut_stat{$base}{$cut}{$mutType} )  { $allMut_stat{$base}{$cut}{$mutType} = $count; }
            else  { $allMut_stat{$base}{$cut}{$mutType} += $count; }
        }
    }

    1;
}

sub outputMutStat 
{
    my $outFile = shift;
    my %parameters = @_;

    my $total = 0; my $percentage = 0; my $printString = "";
    open ( OUT, ">$outFile" );
    print OUT printStatHeader ( "allBase" );
    foreach my $base ( "A", "T", "G", "C" ) {
        ( $printString, $total ) = printAllCount ( $base );
        print OUT $base, $printString;
        $printString = printAllPerc ( $base, $total );
        print OUT $printString;
        if ( defined $parameters{icSHAPE} ) {
            ( $printString, $total ) = printAllCount ( $base, type => "highCut" );
            print OUT $printString;
            $printString = printAllPerc ( $base, $total, type => "highCut" );
            print OUT $printString;
            ( $printString, $total ) = printAllCount ( $base, type => "lowCut" );
            print OUT $printString;
            $printString = printAllPerc ( $base, $total, type => "lowCut" );
            print OUT $printString;
        }
        print OUT "\n";
    }

    foreach my $seqID ( sort {$a cmp $b} ( keys %seqMut_stat ) ) {
        print OUT printStatHeader ( $seqID );
        foreach my $base ( "A", "T", "G", "C" ) {
            ( $printString, $total ) = printSeqCount ( $seqID, $base );
            print OUT $base, $printString;
            $printString = printSeqPerc ( $seqID, $base, $total );
            print OUT $printString;
            if ( defined $parameters{icSHAPE} ) {
                ( $printString, $total ) = printSeqCount ( $seqID, $base, type => "highCut" );
                print OUT $printString;
                $printString = printSeqPerc ( $seqID, $base, $total, type => "highCut" );
                print OUT $printString;
                ( $printString, $total ) = printSeqCount ( $seqID, $base, type => "lowCut" );
                print OUT $printString;
                $printString = printSeqPerc ( $seqID, $base, $total, type => "lowCut" );
                print OUT $printString;
            }
            print OUT "\n";
        }
    }
    close OUT;

    1;
}

sub printStatHeader
{
    my $label = shift;

    my $string = "";
    foreach my $col ( $label, "A", "T", "G", "C", "deletion" ) {  $string .= $col . "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  $string .= $col . "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  $string .= $col . "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  $string .= $col . "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  $string .= $col . "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  $string .= $col . "\t";  }
    $string .= "\n";

    return $string;
}

sub printAllCount
{
    my $base = shift;
    my %parameters = @_;

    my $string = "";
    my $total = 0;
    if ( defined $parameters{type} ) {
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{$parameters{type}}{$mut} ) {
                $total += $allMut_stat{$base}{$parameters{type}}{$mut};
                $string .= "\t" . $allMut_stat{$base}{$parameters{type}}{$mut}; 
            }
            else { $string .="\t0"; }
        }
    }
    else {
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{$mut} ) {
                $total += $allMut_stat{$base}{$mut};
                $string .= "\t" . $allMut_stat{$base}{$mut}; 
            }
            else { $string .="\t0"; }
        }
    }

    return ( $string, $total );
}

sub printSeqCount
{
    my $seqID = shift;
    my $base = shift;
    my %parameters = @_;

    my $string = "";
    my $total = 0;
    if ( defined $parameters{type} ) {
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $seqMut_stat{$seqID}{$base}{$parameters{type}}{$mut} ) {
                $total += $seqMut_stat{$seqID}{$base}{$parameters{type}}{$mut};
                $string .= "\t" . $seqMut_stat{$seqID}{$base}{$parameters{type}}{$mut}; 
            }
            else { $string .="\t0"; }
        }
    }
    else {
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $seqMut_stat{$seqID}{$base}{$mut} ) { 
                $total += $seqMut_stat{$seqID}{$base}{$mut};
                $string .= "\t" . $seqMut_stat{$seqID}{$base}{$mut}; 
            }
            else { $string .="\t0"; }
        }
    }

    return ( $string, $total );
}

sub printAllPerc
{
    my $base = shift; my $total = shift;
    my %parameters = @_;

    my $string = "";
    my $percentage = "-";
    if ( defined $parameters{type} ) {
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{$parameters{type}}{$mut} ) { 
                if ( $base ne $mut )  { $percentage = sprintf ( "%.4f", $allMut_stat{$base}{$parameters{type}}{$mut} / $total ); }
                $string .= "\t". $percentage; 
            }
            else { $string .= "\t0"; }
        }
    }
    else {
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{$mut} ) { 
                if ( $base ne $mut )  { $percentage = sprintf ( "%.4f", $allMut_stat{$base}{$mut} / $total ); }
                $string .= "\t". $percentage; 
            }
            else { $string .= "\t0"; }
        }
    }

    return $string;
}

sub printSeqPerc
{
    my $seqID = shift; my $base = shift; my $total = shift;
    my %parameters = @_;

    my $string = "";
    my $percentage = "-";
    if ( defined $parameters{type} ) {
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $seqMut_stat{$seqID}{$base}{$parameters{type}}{$mut} ) { 
                if ( $base ne $mut )  { $percentage = sprintf ( "%.4f", $seqMut_stat{$seqID}{$base}{$parameters{type}}{$mut} / $total ); }
                $string .= "\t". $percentage; 
            }
            else { $string .= "\t0"; }
        }
    }
    else {
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $seqMut_stat{$seqID}{$base}{$mut} ) { 
                if ( $base ne $mut )  { $percentage = sprintf ( "%.4f", $seqMut_stat{$seqID}{$base}{$mut} / $total ); }
                $string .= "\t". $percentage; 
            }
            else { $string .= "\t0"; }
        }
    }

    return $string;
}



sub readIcSHAPE
{
    my $icSHAPE = shift;
    my %parameters = @_;

    my $count = 0;
    open ( SH, $icSHAPE ) or die ( "Error in reading icSHAPE file $icSHAPE!\n" );
    print STDERR "read icSHAPE file $icSHAPE...\n";
    while ( my $line = <SH> ) {
        $count++;
        chomp $line;

        my ( $id, $length, $rpkm, @scores ) = split ( /\t/, $line );
        $seq_icSHAPE{$id} = \@scores;
    }
    close SH;

    return $count;
}

sub readFasta
{
    my $fasta = shift;
    my %parameters = @_;

    my $count = 0;
    open ( FA, $fasta ) or die ( "Error in reading fasta file $fasta!\n" );
    print STDERR "read fasta file $fasta...\n";
    if ( $parameters{simple} ) {
        while ( my $line = <FA> ) {
            $count++;
            chomp $line;
            my $id = substr ( $line, 1 );
            $id =~ s/^(\s+)//g;

            $line = <FA>;
            chomp $line;
            $seq_fasta{$id} = $line;
        }
    }
    close FA;

    return $count;
}

sub _localReAlignment  
{
    my ( $seqID, $pos, $op ) = @_;

    my $match = 0;
    my $offset = 0;
    my $deletion = uc ( substr ( $op, 1 ) );
    my $substr = uc ( substr ( $seq_fasta{$seqID}, $pos-length($deletion)-1, 3*length($deletion) ) );
    my $result = index($substr, $deletion, $offset);
    while ($result != -1) {
        $match++;
        $offset = $result + 1;
        $result = index($substr, $deletion, $offset);
    }

    print $deletion, "\t", $substr, "\t", $match, "\n" if ( $_debug );
    return $match;
}

sub _append
{
    my $string = shift;
    my $add = shift;

    if ( defined $string ) {  $add = $string . ";" . $add; }
    return $add;
}

sub _parseMD
{
    my $md = shift;
    my %parameters = @_;

    my @op = split ( /[0-9]+/, $md );                   # MD: [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
    my @opSize = split ( /[A-Z^]+/, $md );

    if ( $_debug ) {
        print STDERR "MD\t", $md, "\nopOper";
        _printArray ( \@op );
        print STDERR "\nopSize";
        _printArray ( \@opSize );
        print STDERR "\n";
    }

    return ( \@op, \@opSize );

    1;
}


sub _parseCigar
{
    my $cigar = shift;
    my %parameters = @_;

    my @match = split ( /[0-9]+/, $cigar ); shift @match; # CIGAR: \*|([0-9]+[MIDNSHPX=])+
    my @matchSize = split ( /[MIDNSHPX=]/, $cigar );

    if ( $_debug ) {
        print STDERR "CIGAR\t", $cigar, "\nmatchOper";
        _printArray ( \@match );
        print STDERR "\nmatchSize";
        _printArray ( \@matchSize );
        print STDERR "\n";
    }

    if ( $parameters{getLargestM} ) {
        my $largestM = 0;
        for ( my $idx = 0; $idx < scalar ( @match ); $idx++ ) 
        { if ( $match[$idx] eq "M" ) { if ( $matchSize[$idx] > $largestM ) { $largestM = $matchSize[$idx]; } } }

        return $largestM;
    }
    if ( $parameters{getLeadingS} ) {
        if ( $match[0] ne "S" ) { print STDERR "Warning! unexpected CIGAR string: $cigar!\n"; return 0; }
        else { return $matchSize[0]; }
    }
    elsif ( $parameters{getMatchLen} ) { }
    else { return ( \@match, \@matchSize ); }

    1;
}

sub _printArray
{
    my $ref_array = shift;

    for ( my $idx = 0; $idx < scalar ( @{$ref_array} ); $idx++ ) {
        if ( not defined $ref_array->[$idx] ) { print STDERR "\tnotDefined"; }
        elsif ( $ref_array->[$idx] eq "" ) { print STDERR "\tblank"; }
        else { print STDERR "\t", $ref_array->[$idx]; }
    }

    1;
}

