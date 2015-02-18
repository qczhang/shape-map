#! /usr/bin/perl
#
use strict;
use warnings;

my $_debug = 0;

use Data::Dumper;

my $inputSam = shift;
my $outputMut = shift;
my $refSeq = "/home/qczhang/database/ensembl/current/mouse/gtf/transcriptome.fa";
my $icSHAPE = "/home/qczhang/shape-seq/new/analysis/all.polya/LIB_NAI-LIB_DMSO.PolyA.invitro.valid.enrich";

&main ( $inputSam, $outputMut, icSHAPE => $icSHAPE, refSeq => $refSeq );

sub main
{
    &init ();

    my $samFileList = shift;
    my $mutFile = shift;
    my %parameters = @_;

    my %seq = (); #&readFasta ( $parameters{refSeq}, \%seq, simple => 1 );
    my %seq_icSHAPE = (); #&readIcSHAPE ( $parameters{icSHAPE}, \%seq_icSHAPE );

    my %seq_mut = ();
    my @samFiles = split ( /:/, $samFileList );
    foreach my $samFile ( @samFiles ) { &readSam ( $samFile, \%seq, \%seq_mut ); }
    print Dumper \%seq_mut;
    exit;

    my %seqMut_stat = (); my %allMut_stat = ();
    &mutStatistics ( \%seq, \%seq_icSHAPE, \%seq_mut, \%seqMut_stat, \%allMut_stat, lowCut => 0.05, highCut => 0.4 );

    &outputMutStat ( \%seqMut_stat, \%allMut_stat, $mutFile );
}

    if ( not $inputSam ) { die "Usage: $0 input_sam output_stat fasta shape\n"; }

sub init
{
    if ( not $inputSam ) { die "Usage: $0 input_sam output_stat fasta shape\n"; }
    if ( ( not defined $outputMut ) or ( $outputMut eq "NULL" ) ) {  $outputMut = "output.stat";  }
    if ( ( not defined $refSeq ) or ( $refSeq eq "NULL" ) ) { $refSeq = "/home/qczhang/database/ensembl/current/mouse/gtf/transcriptome.fa"; }
    if ( ( not defined $icSHAPE ) or ( $refSeq eq "NULL" ) ) { $icSHAPE = "/home/qczhang/shape-seq/new/analysis/all.polya/LIB_NAI-LIB_DMSO.PolyA.invitro.valid.enrich"; }
}

sub readSam
{
    my $samFile = shift;
    my $ref_seq = shift;
    my $ref_seq_mut = shift;

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
                my $mdString = ""; for ( my $idx = 11; $idx < scalar ( @data ); $idx++ ) { if ( $data[$idx] =~ /MD:/ ) { $mdString = substr ( $data[$idx], 5 ); last; } }
                &parseMut ( $ref_seq, $ref_seq_mut, $data[2], $data[3], $data[5], $data[9],$mdString );
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
                &parseMut2 ( $ref_seq_mut, $data[2], $data[3], $data2[3], $data[5], $data2[5], $data[9], $data2[9], $mdString, $mdString2 );
            }
        }
    }
    close SAM;

    return $lineCount;
}

sub parseMut
{
    my $ref_seq = shift; my $ref_mut = shift;
    my $seqID = shift; my $pos = shift; my $cigar = shift; my $readSeq = shift; my $md = shift;

    my ( $ref_match, $ref_matchSize ) = _parseCigar ( $cigar );
    my ( $ref_op, $ref_opSize ) = _parseMD ( $md );

    my @op = split ( /[0-9]+/, $md );                   # MD: [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
    my @opsize = split ( /[A-Z^]+/, $md );
    if ( $_debug ) {
        print "CIGAR\t", $cigar, "\nmatch";
        _printString ( \@match );
        print "\nmatchSize";
        _printString ( \@matchSize );
        print "\n";
        print "MD\t", $md, "\nop";
        _printString ( \@op );
        print "\nopsize";
        _printString ( \@opsize );
        print "\n";
    }

    my $refPos = $pos - 1;
    my $readPos = 0;  
    my $removeClip = 0; my @insertions = (); my @insertionSize = ();
    for ( my $idx = 1; $idx < scalar ( @match ); $idx++ ) {
        ## cigar string are always numbers-matchsymbols, so $match[1] corresponds to $matchSize[0]
        if ( ( $match[$idx] eq "S" ) or ( $match[$idx] eq "H" ) ) { 
            if ( $idx == 1 ) { 
#                $ref_mut->{$seqID}{"clip"}[$refPos] = _append ( $ref_mut->{$seqID}{"clip"}[$refPos], substr ( $readSeq, $readPos, $matchSize[$idx-1] ) . "-" );
#                if ( $match[$idx] eq "H" ) { $ref_mut->{$seqID}{"clip"}[$refPos-1] = _append ( $ref_mut->{$seqID}{"clip"}[$refPos], "csna" ); }
                $removeClip = $matchSize[$idx-1];
            }
            elsif ( $idx == $#match ) { 
#                $ref_mut->{$seqID}{"clip"}[$refPos] = _append ( $ref_mut->{$seqID}{"clip"}[$refPos], "-" . substr ( $readSeq, $readPos, $matchSize[$idx-1] ) );
#                if ( $match[$idx] eq "H" ) { $ref_mut->{$seqID}{"clip"}[$refPos] = _append ( $ref_mut->{$seqID}{"clip"}[$refPos], "csna" ); }
            }
            else {
                print STDERR "Warning! clipping inside an alignment!\n";
                print STDERR "\t$cigar\t$readSeq\n";
            }

            $readPos += $matchSize[$idx-1];
        }
        elsif ( $match[$idx] eq "I" ) {
#            $ref_mut->{$seqID}{"insertion"}[$refPos] = _append ( $ref_mut->{$seqID}{"insertion"}[$refPos], substr ( $readSeq, $readPos, $matchSize[$idx-1] ) );
            $readPos += $matchSize[$idx-1];
            push ( @insertions, $refPos+1 );
            push ( @insertionSize, $matchSize[$idx-1] );
        }
        elsif ( $match[$idx] eq "P" ) {
        }
        elsif ( $match[$idx] eq "D" ) {
            $refPos += $matchSize[$idx-1];     
            # $ref_mut->{$seqID}{"deletion"}[$refPos-1] = _append ( $ref_mut->{$seqID}{"deletion"}[$refPos-1], $matchSize[$idx-1] );
        }
        elsif ( $match[$idx] eq "X" ) {
            ## not tested
            #$ref_mut->{$seqID}{"mismatch"}[$refPos] = _append ( $ref_mut->{$seqID}{"mismatch"}[$refPos], substr ( $readSeq, $readPos, $matchSize[$idx-1] ) );
            $refPos += $matchSize[$idx-1];     
            $readPos += $matchSize[$idx-1];
        }
        elsif ( ( $match[$idx] eq "=" ) or ( $match[$idx] eq "M" ) ) {
            $refPos += $matchSize[$idx-1];     
            $readPos += $matchSize[$idx-1];
        }
    }

    $refPos = $pos;                 ## point to 1-base index of next event
    $readPos = $removeClip + 1;     ## the same

    if ( $op[0] ) {
        if ( $op[0] =~ /^\^/ ) {
            ## removal of ambiguity
            my $delLen = length ( $op[0] ) - 1;
            for ( my $idxPos = 0; $idxPos < $delLen-1; $idxPos++ ) { $ref_mut->{$seqID}{"same"}[$refPos+$idxPos-1]++; }
            $refPos += $delLen;
            if ( ( defined $op[1] ) and ( not defined $opsize[1] ) ) { print STDERR "error!\n"; return -1; }

            if ( ( defined $op[1] ) and ( $op[1] ) and ( not $opsize[0] ) ) { $ref_mut->{$seqID}{"same"}[$refPos-1-1]++; }
            else {
                my $matchInLocal = _localReAlignment ( $ref_seq, $seqID, $refPos-$delLen, $op[0] );
                print "multiple matching? - $matchInLocal\n" if ( $_debug );
                if ( $matchInLocal == 1 ) {
                    print "\tmatch only once, no realignment needed.\n" if ( $_debug );
                    $ref_mut->{$seqID}{"deletion"}[$refPos-1-1] = _append ( $ref_mut->{$seqID}{"deletion"}[$refPos-1-1], $op[0] );
                    ## you want to label the deletion even at the last base of deleted fragment
                }
                else { $ref_mut->{$seqID}{"same"}[$refPos-1-1]++; }
            }
        }
        elsif ( $op[0] =~ /[A-Z]/ ) {
            $refPos++; $readPos++;
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            if ( ( defined $op[1] ) and ( not defined $opsize[1] ) ) { print STDERR "error!\n"; return -1; }

            if ( ( defined $op[1] ) and ( $op[1] ) and ( not $opsize[1] ) ) { $ref_mut->{$seqID}{"same"}[$refPos-1-1]++; }
            else { $ref_mut->{$seqID}{"mutation"}[$refPos-1-1] = _append ( $ref_mut->{$seqID}{"mutation"}[$refPos-1-1], substr ( $readSeq, $readPos-1-1, 1 ) ); }
        }

        $op[0] = "";
        shift @opsize;
    }

    for ( my $idx = 0; $idx <= scalar ( @op ); $idx++ ) {
        ## make sure the order is not affected by the real operation strings!!!!!!
        if ( ( not defined $op[$idx] ) or ( not $op[$idx] ) ) {
        }
        elsif ( $op[$idx] =~ /^\^/ ) {
            ## removal of ambiguity
            my $delLen = length ( $op[$idx] ) - 1;
            for ( my $idxPos = 0; $idxPos < $delLen-1; $idxPos++ ) { $ref_mut->{$seqID}{"same"}[$refPos+$idxPos-1]++; }
            $refPos += $delLen;
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            if ( ( defined $op[$idx+1] ) and ( $op[$idx+1] ) and ( not $opsize[$idx] ) ) { $ref_mut->{$seqID}{"same"}[$refPos-1-1]++; }
            else {
                my $matchInLocal = _localReAlignment ( $ref_seq, $seqID, $refPos-$delLen, $op[$idx] );
                print "multiple matching? - $matchInLocal\n" if ( $_debug );
                if ( $matchInLocal == 1 ) {
                    print "\tmatch only once, no realignment needed.\n" if ( $_debug );
                    $ref_mut->{$seqID}{"deletion"}[$refPos-1-1] = _append ( $ref_mut->{$seqID}{"deletion"}[$refPos-1-1], $op[$idx] );
                    ## you want to label the deletion even at the last base of deleted fragment
                }
                else { $ref_mut->{$seqID}{"same"}[$refPos-1-1]++; }
            }
        }
        elsif ( $op[$idx] =~ /[A-Z]/ ) {
            $refPos++; $readPos++;
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            if ( ( defined $op[$idx+1] ) and ( $op[$idx+1] ) and ( not $opsize[$idx] ) ) { $ref_mut->{$seqID}{"same"}[$refPos-1-1]++; }
            else { $ref_mut->{$seqID}{"mutation"}[$refPos-1-1] = _append ( $ref_mut->{$seqID}{"mutation"}[$refPos-1-1], substr ( $readSeq, $readPos-1-1, 1 ) ); }
        }

        if ( $opsize[$idx] ) {
            for ( my $idxPos = 0; $idxPos < $opsize[$idx]; $idxPos++ ) { $ref_mut->{$seqID}{"same"}[$refPos+$idxPos-1]++; }

            $refPos += $opsize[$idx];
            $readPos += $opsize[$idx];

            while ( ( scalar (@insertions ) ) and ( $refPos >= $insertions[0] ) )  {
                $readPos += $insertionSize[0];
                shift ( @insertions );
                shift ( @insertionSize );
            }
        }
    }


    1;
}

sub _parseMut2
{
    my $ref_mut = shift;

    my $seqID = shift;
    my $pos = shift;
    my $pos2 = shift;
    my $cigar = shift;
    my $cigar2 = shift;
    my $readSeq = shift;
    my $readSeq2 = shift;
    my $md = shift;
    my $md2 = shift;

    my @match = split ( /[0-9]+/, $cigar );             # CIGAR: \*|([0-9]+[MIDNSHPX=])+
    my @matchSize = split ( /[MIDNSHPX=]/, $cigar );
    my @op = split ( /[0-9]+/, $md );                   # MD: [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
    my @opsize = split ( /[A-Z^]+/, $md );
    my @match2 = split ( /[0-9]+/, $cigar2 );             # CIGAR: \*|([0-9]+[MIDNSHPX=])+
    my @matchSize2 = split ( /[MIDNSHPX=]/, $cigar2 );
    my @op2 = split ( /[0-9]+/, $md2 );                   # MD: [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
    my @opsize2 = split ( /[A-Z^]+/, $md2 );

    my $refPos = $pos - 1;
    my $readPos = 0;  
    my $removeClip = 0; my @insertions = (); my @insertionSize = ();
    for ( my $idx = 1; $idx < scalar ( @match ); $idx++ ) {
        ## cigar string are always numbers-matchsymbols, so $match[1] corresponds to $matchSize[0]
        if ( ( $match[$idx] eq "S" ) or ( $match[$idx] eq "H" ) ) { 
            if ( $idx == 1 ) { 
#                $ref_mut->{$seqID}{"clip"}[$refPos] = _append ( $ref_mut->{$seqID}{"clip"}[$refPos], substr ( $readSeq, $readPos, $matchSize[$idx-1] ) . "-" );
#                if ( $match[$idx] eq "H" ) { $ref_mut->{$seqID}{"clip"}[$refPos-1] = _append ( $ref_mut->{$seqID}{"clip"}[$refPos], "csna" ); }
                $removeClip = $matchSize[$idx-1];
            }
            elsif ( $idx == $#match ) { 
#                $ref_mut->{$seqID}{"clip"}[$refPos] = _append ( $ref_mut->{$seqID}{"clip"}[$refPos], "-" . substr ( $readSeq, $readPos, $matchSize[$idx-1] ) );
#                if ( $match[$idx] eq "H" ) { $ref_mut->{$seqID}{"clip"}[$refPos] = _append ( $ref_mut->{$seqID}{"clip"}[$refPos], "csna" ); }
            }
            else {
                print STDERR "Warning! clipping inside an alignment!\n";
                print STDERR "\t$cigar\t$readSeq\n";
            }

            $readPos += $matchSize[$idx-1];
        }
        elsif ( $match[$idx] eq "I" ) {
#            $ref_mut->{$seqID}{"insertion"}[$refPos] = _append ( $ref_mut->{$seqID}{"insertion"}[$refPos], substr ( $readSeq, $readPos, $matchSize[$idx-1] ) );
            $readPos += $matchSize[$idx-1];
            push ( @insertions, $refPos+1 );
            push ( @insertionSize, $matchSize[$idx-1] );
        }
        elsif ( $match[$idx] eq "P" ) {
        }
        elsif ( $match[$idx] eq "D" ) {
            $refPos += $matchSize[$idx-1];     
            # $ref_mut->{$seqID}{"deletion"}[$refPos-1] = _append ( $ref_mut->{$seqID}{"deletion"}[$refPos-1], $matchSize[$idx-1] );
        }
        elsif ( $match[$idx] eq "X" ) {
            ## not tested
            #$ref_mut->{$seqID}{"mismatch"}[$refPos] = _append ( $ref_mut->{$seqID}{"mismatch"}[$refPos], substr ( $readSeq, $readPos, $matchSize[$idx-1] ) );
            $refPos += $matchSize[$idx-1];     
            $readPos += $matchSize[$idx-1];
        }
        elsif ( ( $match[$idx] eq "=" ) or ( $match[$idx] eq "M" ) ) {
            $refPos += $matchSize[$idx-1];     
            $readPos += $matchSize[$idx-1];
        }
    }

    $refPos = $pos;                 ## point to 1-base index of next event
    $readPos = $removeClip + 1;     ## the same
    for ( my $idx = 0; $idx < scalar ( @op ); $idx++ ) {
        if ( not $op[$idx] ) {
        }
        elsif ( $op[$idx] =~ /^\^/ ) {
            my $delLen = length ( $op[$idx] ) - 1;
            $refPos += $delLen;
            $ref_mut->{$seqID}{"deletion"}[$refPos-1-1] = _append ( $ref_mut->{$seqID}{"deletion"}[$refPos-1-1], $op[$idx] );
            ## you want to label the deletion even at the last base of deleted fragment
        }
        elsif ( $op[$idx] =~ /[A-Z]/ ) {
            $refPos++;
            $readPos++;
            $ref_mut->{$seqID}{"mutation"}[$refPos-1-1] = _append ( $ref_mut->{$seqID}{"mutation"}[$refPos-1-1], substr ( $readSeq, $readPos-1-1, 1 ) );
        }

        if ( defined $opsize[$idx] ) {
            $refPos += $opsize[$idx];
            $readPos += $opsize[$idx];

            while ( ( scalar (@insertions ) ) and ( $refPos >= $insertions[0] ) )  {
                $readPos += $insertionSize[0];
                shift ( @insertions );
                shift ( @insertionSize );
            }
        }
    }

    1;
}


sub mutStatistics
{
    my $ref_seq = shift;
    my $ref_seqIcSHAPE = shift;
    my $ref_mut = shift;
    my $ref_seqMutStat = shift;
    my $ref_allMutStat = shift;
    my %parameters = @_;

    foreach my $seqID ( keys %{$ref_mut} ) {
        if ( not defined $ref_seq->{$seqID} ) {
            print STDERR "ERROR! Sequence of $seqID not found. ...skipped\n";
            next;
        }

        my $seq = uc ( $ref_seq->{$seqID} );
        for ( my $idx = 0; $idx < length ( $seq ); $idx++ ) {
            my $base = substr ( $seq, $idx, 1 );
            my $icSHAPE = ( defined $ref_seqIcSHAPE->{$seqID} ) ? $ref_seqIcSHAPE->{$seqID}[$idx] : "NULL"; 
            if ( defined $ref_mut->{$seqID}{same}[$idx] ) {
                if ( not defined $ref_seqMutStat->{$seqID}{$base}{$base} )  { $ref_seqMutStat->{$seqID}{$base}{$base} = $ref_mut->{$seqID}{same}[$idx]; }
                else  { $ref_seqMutStat->{$seqID}{$base}{$base} += $ref_mut->{$seqID}{same}[$idx]; }
                if ( not defined $ref_allMutStat->{$base}{$base} )  { $ref_allMutStat->{$base}{$base} = $ref_mut->{$seqID}{same}[$idx]; }
                else  { $ref_allMutStat->{$base}{$base} += $ref_mut->{$seqID}{same}[$idx]; }

                if ( $icSHAPE ne "NULL" ) {
                    if ( $icSHAPE < $parameters{lowCut} ) {
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{lowCut}{$base} )  { $ref_seqMutStat->{$seqID}{$base}{lowCut}{$base} = $ref_mut->{$seqID}{same}[$idx]; }
                        else  { $ref_seqMutStat->{$seqID}{$base}{lowCut}{$base} += $ref_mut->{$seqID}{same}[$idx]; }
                        if ( not defined $ref_allMutStat->{$base}{lowCut}{$base} )  { $ref_allMutStat->{$base}{lowCut}{$base} = $ref_mut->{$seqID}{same}[$idx]; }
                        else  { $ref_allMutStat->{$base}{lowCut}{$base} += $ref_mut->{$seqID}{same}[$idx]; }
                    }
                    elsif ( $icSHAPE > $parameters{highCut} ) {
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{highCut}{$base} )  { $ref_seqMutStat->{$seqID}{$base}{highCut}{$base} = $ref_mut->{$seqID}{same}[$idx]; }
                        else  { $ref_seqMutStat->{$seqID}{$base}{highCut}{$base} += $ref_mut->{$seqID}{same}[$idx]; }
                        if ( not defined $ref_allMutStat->{$base}{highCut}{$base} )  { $ref_allMutStat->{$base}{highCut}{$base} = $ref_mut->{$seqID}{same}[$idx]; }
                        else  { $ref_allMutStat->{$base}{highCut}{$base} += $ref_mut->{$seqID}{same}[$idx]; }
                    }
                }
            }
            if ( defined $ref_mut->{$seqID}{deletion}[$idx] ) {
                my $insD = () = $ref_mut->{$seqID}{deletion}[$idx] =~ /;/gi; $insD++;

                if ( not defined $ref_seqMutStat->{$seqID}{$base}{deletion} )  { $ref_seqMutStat->{$seqID}{$base}{deletion} = $insD; }
                else { $ref_seqMutStat->{$seqID}{$base}{deletion} += $insD; }
                if ( not defined $ref_allMutStat->{$base}{deletion} )  { $ref_allMutStat->{$base}{deletion} = $insD; }
                else { $ref_allMutStat->{$base}{deletion} += $insD; }
                if ( $icSHAPE ne "NULL" ) {
                    if ( $icSHAPE < $parameters{lowCut} ) {
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{lowCut}{deletion} )  { $ref_seqMutStat->{$seqID}{$base}{lowCut}{deletion} = $insD; }
                        else { $ref_seqMutStat->{$seqID}{$base}{lowCut}{deletion} += $insD; }
                        if ( not defined $ref_allMutStat->{$base}{lowCut}{deletion} )  { $ref_allMutStat->{$base}{lowCut}{deletion} = $insD; }
                        else { $ref_allMutStat->{$base}{lowCut}{deletion} += $insD; }
                    }
                    elsif ( $icSHAPE > $parameters{highCut} ) {
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{highCut}{deletion} )  { $ref_seqMutStat->{$seqID}{$base}{highCut}{deletion} = $insD; }
                        else { $ref_seqMutStat->{$seqID}{$base}{highCut}{deletion} += $insD; }
                        if ( not defined $ref_allMutStat->{$base}{highCut}{deletion} )  { $ref_allMutStat->{$base}{highCut}{deletion} = $insD; }
                        else { $ref_allMutStat->{$base}{highCut}{deletion} += $insD; }
                    }
                }
            }
            if ( defined $ref_mut->{$seqID}{mutation}[$idx] ) {
                my $insA = () = $ref_mut->{$seqID}{mutation}[$idx] =~ /A/gi;
                my $insT = () = $ref_mut->{$seqID}{mutation}[$idx] =~ /T/gi;
                my $insG = () = $ref_mut->{$seqID}{mutation}[$idx] =~ /G/gi;
                my $insC = () = $ref_mut->{$seqID}{mutation}[$idx] =~ /C/gi;

                if ( not defined $ref_seqMutStat->{$seqID}{$base}{A} )  { $ref_seqMutStat->{$seqID}{$base}{A} = $insA; }
                else { $ref_seqMutStat->{$seqID}{$base}{A} += $insA; }
                if ( not defined $ref_seqMutStat->{$seqID}{$base}{T} )  { $ref_seqMutStat->{$seqID}{$base}{T} = $insT; }
                else { $ref_seqMutStat->{$seqID}{$base}{T} += $insT; }
                if ( not defined $ref_seqMutStat->{$seqID}{$base}{G} )  { $ref_seqMutStat->{$seqID}{$base}{G} = $insG; }
                else { $ref_seqMutStat->{$seqID}{$base}{G} += $insG; }
                if ( not defined $ref_seqMutStat->{$seqID}{$base}{C} )  { $ref_seqMutStat->{$seqID}{$base}{C} = $insC; }
                else { $ref_seqMutStat->{$seqID}{$base}{C} += $insC; }

                if ( not defined $ref_allMutStat->{$base}{A} )  { $ref_allMutStat->{$base}{A} = $insA; }
                else { $ref_allMutStat->{$base}{A} += $insA; }
                if ( not defined $ref_allMutStat->{$base}{T} )  { $ref_allMutStat->{$base}{T} = $insT; }
                else { $ref_allMutStat->{$base}{T} += $insT; }
                if ( not defined $ref_allMutStat->{$base}{G} )  { $ref_allMutStat->{$base}{G} = $insG; }
                else { $ref_allMutStat->{$base}{G} += $insG; }
                if ( not defined $ref_allMutStat->{$base}{C} )  { $ref_allMutStat->{$base}{C} = $insC; }
                else { $ref_allMutStat->{$base}{C} += $insC; }

                if ( $icSHAPE ne "NULL" ) {
                    if ( $icSHAPE < $parameters{lowCut} ) {
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{lowCut}{A} )  { $ref_seqMutStat->{$seqID}{$base}{lowCut}{A} = $insA; }
                        else { $ref_seqMutStat->{$seqID}{$base}{lowCut}{A} += $insA; }
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{lowCut}{T} )  { $ref_seqMutStat->{$seqID}{$base}{lowCut}{T} = $insT; }
                        else { $ref_seqMutStat->{$seqID}{$base}{lowCut}{T} += $insT; }
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{lowCut}{G} )  { $ref_seqMutStat->{$seqID}{$base}{lowCut}{G} = $insG; }
                        else { $ref_seqMutStat->{$seqID}{$base}{lowCut}{G} += $insG; }
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{lowCut}{C} )  { $ref_seqMutStat->{$seqID}{$base}{lowCut}{C} = $insC; }
                        else { $ref_seqMutStat->{$seqID}{$base}{lowCut}{C} += $insC; }

                        if ( not defined $ref_allMutStat->{$base}{lowCut}{A} )  { $ref_allMutStat->{$base}{lowCut}{A} = $insA; }
                        else { $ref_allMutStat->{$base}{lowCut}{A} += $insA; }
                        if ( not defined $ref_allMutStat->{$base}{lowCut}{T} )  { $ref_allMutStat->{$base}{lowCut}{T} = $insT; }
                        else { $ref_allMutStat->{$base}{lowCut}{T} += $insT; }
                        if ( not defined $ref_allMutStat->{$base}{lowCut}{G} )  { $ref_allMutStat->{$base}{lowCut}{G} = $insG; }
                        else { $ref_allMutStat->{$base}{lowCut}{G} += $insG; }
                        if ( not defined $ref_allMutStat->{$base}{lowCut}{C} )  { $ref_allMutStat->{$base}{lowCut}{C} = $insC; }
                        else { $ref_allMutStat->{$base}{lowCut}{C} += $insC; }
                    }
                    elsif ( $icSHAPE > $parameters{highCut} ) {
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{highCut}{A} )  { $ref_seqMutStat->{$seqID}{$base}{highCut}{A} = $insA; }
                        else { $ref_seqMutStat->{$seqID}{$base}{highCut}{A} += $insA; }
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{highCut}{T} )  { $ref_seqMutStat->{$seqID}{$base}{highCut}{T} = $insT; }
                        else { $ref_seqMutStat->{$seqID}{$base}{highCut}{T} += $insT; }
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{highCut}{G} )  { $ref_seqMutStat->{$seqID}{$base}{highCut}{G} = $insG; }
                        else { $ref_seqMutStat->{$seqID}{$base}{highCut}{G} += $insG; }
                        if ( not defined $ref_seqMutStat->{$seqID}{$base}{highCut}{C} )  { $ref_seqMutStat->{$seqID}{$base}{highCut}{C} = $insC; }
                        else { $ref_seqMutStat->{$seqID}{$base}{highCut}{C} += $insC; }

                        if ( not defined $ref_allMutStat->{$base}{highCut}{A} )  { $ref_allMutStat->{$base}{highCut}{A} = $insA; }
                        else { $ref_allMutStat->{$base}{highCut}{A} += $insA; }
                        if ( not defined $ref_allMutStat->{$base}{highCut}{T} )  { $ref_allMutStat->{$base}{highCut}{T} = $insT; }
                        else { $ref_allMutStat->{$base}{highCut}{T} += $insT; }
                        if ( not defined $ref_allMutStat->{$base}{highCut}{G} )  { $ref_allMutStat->{$base}{highCut}{G} = $insG; }
                        else { $ref_allMutStat->{$base}{highCut}{G} += $insG; }
                        if ( not defined $ref_allMutStat->{$base}{highCut}{C} )  { $ref_allMutStat->{$base}{highCut}{C} = $insC; }
                        else { $ref_allMutStat->{$base}{highCut}{C} += $insC; }
                    }
                }
            }
        }
    }

    1;
}

sub outputMutStat 
{
    my $ref_seqMutStat = shift;
    my $ref_allMutStat = shift;
    my $outFile = shift;
    my %parameters = @_;

    my $total = 0; my $percentage = 0;
    open ( OUT, ">$outFile" );
    foreach my $col ( "allBase", "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
    foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
    print OUT "\n";
    foreach my $base ( "A", "T", "G", "C" ) {
        print OUT $base;
        $total = 0;
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $ref_allMutStat->{$base}{$mut} ) {
                $total += $ref_allMutStat->{$base}{$mut};
                print OUT "\t", $ref_allMutStat->{$base}{$mut}; 
            }
            else { print OUT "\t0"; }
        }
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $ref_allMutStat->{$base}{$mut} ) { 
                if ( $base eq $mut ) { $percentage = "-" ; }
                else { $percentage = sprintf ( "%.5f", $ref_allMutStat->{$base}{$mut} / $total ); }
                print OUT "\t", $percentage; 
            }
            else { print OUT "\t0"; }
        }

        $total = 0;
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $ref_allMutStat->{$base}{lowCut}{$mut} ) {
                $total += $ref_allMutStat->{$base}{lowCut}{$mut};
                print OUT "\t", $ref_allMutStat->{$base}{lowCut}{$mut}; 
            }
            else { print OUT "\t0"; }
        }
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $ref_allMutStat->{$base}{lowCut}{$mut} ) { 
                if ( $base eq $mut ) { $percentage = "-" ; }
                else { $percentage = sprintf ( "%.5f", $ref_allMutStat->{$base}{lowCut}{$mut} / $total ); }
                print OUT "\t", $percentage; 
            }
            else { print OUT "\t0"; }
        }

        $total = 0;
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $ref_allMutStat->{$base}{highCut}{$mut} ) {
                $total += $ref_allMutStat->{$base}{highCut}{$mut};
                print OUT "\t", $ref_allMutStat->{$base}{highCut}{$mut}; 
            }
            else { print OUT "\t0"; }
        }
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $ref_allMutStat->{$base}{highCut}{$mut} ) { 
                if ( $base eq $mut ) { $percentage = "-" ; }
                else { $percentage = sprintf ( "%.5f", $ref_allMutStat->{$base}{highCut}{$mut} / $total ); }
                print OUT "\t", $percentage; 
            }
            else { print OUT "\t0"; }
        }
        print OUT "\n";
    }

    foreach my $seqID ( sort {$a cmp $b} ( keys %{$ref_seqMutStat} ) ) {
        foreach my $col ( $seqID, "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
        foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
        foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
        foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
        foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
        foreach my $col ( "A", "T", "G", "C", "deletion" ) {  print OUT $col, "\t";  }
        print OUT "\n";
        foreach my $base ( "A", "T", "G", "C" ) {
            print OUT $base;
            $total = 0;
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $ref_seqMutStat->{$seqID}{$base}{$mut} ) { 
                    $total += $ref_seqMutStat->{$seqID}{$base}{$mut};
                    print OUT "\t", $ref_seqMutStat->{$seqID}{$base}{$mut}; 
                }
                else { print OUT "\t0"; }
            }
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $ref_seqMutStat->{$seqID}{$base}{$mut} ) { 
                    if ( $base eq $mut ) { $percentage = "-" ; }
                    else { $percentage = sprintf ( "%.5f", $ref_seqMutStat->{$seqID}{$base}{$mut} / $total ); }
                    print OUT "\t", $percentage; 
                }
                else { print OUT "\t0"; }
            }

            $total = 0;
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $ref_seqMutStat->{$seqID}{$base}{lowCut}{$mut} ) { 
                    $total += $ref_seqMutStat->{$seqID}{$base}{lowCut}{$mut};
                    print OUT "\t", $ref_seqMutStat->{$seqID}{$base}{lowCut}{$mut}; 
                }
                else { print OUT "\t0"; }
            }
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $ref_seqMutStat->{$seqID}{$base}{lowCut}{$mut} ) { 
                    if ( $base eq $mut ) { $percentage = "-" ; }
                    else { $percentage = sprintf ( "%.5f", $ref_seqMutStat->{$seqID}{$base}{lowCut}{$mut} / $total ); }
                    print OUT "\t", $percentage; 
                }
                else { print OUT "\t0"; }
            }

            $total = 0;
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $ref_seqMutStat->{$seqID}{$base}{highCut}{$mut} ) { 
                    $total += $ref_seqMutStat->{$seqID}{$base}{highCut}{$mut};
                    print OUT "\t", $ref_seqMutStat->{$seqID}{$base}{highCut}{$mut}; 
                }
                else { print OUT "\t0"; }
            }
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $ref_seqMutStat->{$seqID}{$base}{highCut}{$mut} ) { 
                    if ( $base eq $mut ) { $percentage = "-" ; }
                    else { $percentage = sprintf ( "%.5f", $ref_seqMutStat->{$seqID}{$base}{highCut}{$mut} / $total ); }
                    print OUT "\t", $percentage; 
                }
                else { print OUT "\t0"; }
            }
            print OUT "\n";
        }
    }

    1;
}

sub readIcSHAPE
{
    my $icSHAPE = shift;
    my $ref_icSHAPE = shift;
    my %parameters = @_;

    my $count = 0;
    open ( SH, $icSHAPE ) or die ( "Error in reading icSHAPE file $icSHAPE!\n" );
    print STDERR "read icSHAPE file $icSHAPE...\n";
    while ( my $line = <SH> ) {
        $count++;
        chomp $line;

        my ( $id, $length, $rpkm, @scores ) = split ( /\t/, $line );
        $ref_icSHAPE->{$id} = \@scores;
    }
    close SH;

    return $count;
}

sub readFasta
{
    my $fasta = shift;
    my $ref_seq = shift;
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
            $ref_seq->{$id} = $line;
        }
    }
    close FA;

    return $count;
}

sub _localReAlignment  
{
    my ( $ref_seq, $seqID, $pos, $op ) = @_;

    my $match = 0;
    my $offset = 0;
    my $deletion = uc ( substr ( $op, 1 ) );
    my $substr = uc ( substr ( $ref_seq->{$seqID}, $pos-length($deletion)-1, 3*length($deletion) ) );
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

sub _printString
{
    my $ref_array = shift;

    for ( my $idx = 0; $idx < scalar ( @{$ref_array} ); $idx++ ) {
        if ( not defined $ref_array->[$idx] ) { print "\tn.d."; }
        elsif ( $ref_array->[$idx] eq "" ) { print "\tblank"; }
        else { print "\t", $ref_array->[$idx]; }
    }

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
        for ( my $idx = 0; $idx < scalar ( @match ); $idx++ ) {
            if ( $match[$idx] eq "M" ) {
                if ( $matchSize[$idx] > $largestM ) {
                    $largestM = $matchSize[$idx];
                }
            }
        }

        return $largestM;
    }
    if ( $parameters{getLeadingS} ) {
        if ( $match[0] ne "S" ) {
            print STDERR "Warning! unexpected CIGAR string: $cigar!\n";
            return 0;
        }
        else { return $matchSize[0]; }
    }
    elsif ( $parameters{getMatchLen} ) {
    }
    else {
        return ( \@match, \@matchSize );
    }

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
        for ( my $idx = 0; $idx < scalar ( @match ); $idx++ ) {
            if ( $match[$idx] eq "M" ) {
                if ( $matchSize[$idx] > $largestM ) {
                    $largestM = $matchSize[$idx];
                }
            }
        }

        return $largestM;
    }
    if ( $parameters{getLeadingS} ) {
        if ( $match[0] ne "S" ) {
            print STDERR "Warning! unexpected CIGAR string: $cigar!\n";
            return 0;
        }
        else { return $matchSize[0]; }
    }
    elsif ( $parameters{getMatchLen} ) {
    }
    else {
        return ( \@match, \@matchSize );
    }

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

