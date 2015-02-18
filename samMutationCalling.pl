#! /usr/bin/perl
#
use strict;
use warnings;

my $_debug = 0;

#use Data::Dumper;

my $inputSam = shift;
my $outputMut = shift;
my $refSeq = "/home/qczhang/database/ensembl/current/mouse/gtf/transcriptome.fa";
my $icSHAPE = "/home/qczhang/shape-seq/new/analysis/all.polya/LIB_NAI-LIB_DMSO.PolyA.invitro.valid.enrich";

my %seq = (); 
my %seq_icSHAPE = ();
my %seq_mut = ();
my %seqMut_stat = (); 
my %allMut_stat = ();

&main ( $inputSam, $outputMut, icSHAPE => $icSHAPE, refSeq => $refSeq );

sub main
{
    &init ();

    my $samFileList = shift;
    my $mutFile = shift;
    my %parameters = @_;

    &readFasta ( $parameters{refSeq}, simple => 1 );
    &readIcSHAPE ( $parameters{icSHAPE}, );

    my @samFiles = split ( /:/, $samFileList );
    foreach my $samFile ( @samFiles ) { &readSam ( $samFile ); }
    print Dumper \%seq_mut;
    exit;

    &mutStatistics ( lowCut => 0.05, highCut => 0.4 );

    &outputMutStat ( $mutFile );
}

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

    my ( $headSoftClip, $ref_insertions, $ref_insertionSize ) = collectAlignInfoCIGAR ( );
    collectAlignInfoMD ( );

    1;
}

sub _parseMut2
{
}

sub mutStatistics
{
    my %parameters = @_;

    foreach my $seqID ( keys %seq_mut ) {
        if ( not defined $seq{$seqID} ) {
            print STDERR "ERROR! Sequence of $seqID not found. ...skipped\n";
            next;
        }

        my $seq = uc ( $seq{$seqID} );
        for ( my $idx = 0; $idx < length ( $seq ); $idx++ ) {
            my $base = substr ( $seq, $idx, 1 );
            my $icSHAPE = ( defined $seq_icSHAPE{$seqID} ) ? $seq_icSHAPE{$seqID}[$idx] : "NULL"; 
            if ( defined $seq_mut{$seqID}{same}[$idx] ) {
                if ( not defined $seqMut_stat{$seqID}{$base}{$base} )  { $seqMut_stat{$seqID}{$base}{$base} = $seq_mut{$seqID}{same}[$idx]; }
                else  { $seqMut_stat{$seqID}{$base}{$base} += $seq_mut{$seqID}{same}[$idx]; }
                if ( not defined $allMut_stat{$base}{$base} )  { $allMut_stat{$base}{$base} = $seq_mut{$seqID}{same}[$idx]; }
                else  { $allMut_stat{$base}{$base} += $seq_mut{$seqID}{same}[$idx]; }

                if ( $icSHAPE ne "NULL" ) {
                    if ( $icSHAPE < $parameters{lowCut} ) {
                        if ( not defined $seqMut_stat{$seqID}{$base}{lowCut}{$base} )  { $seqMut_stat{$seqID}{$base}{lowCut}{$base} = $seq_mut{$seqID}{same}[$idx]; }
                        else  { $seqMut_stat{$seqID}{$base}{lowCut}{$base} += $seq_mut{$seqID}{same}[$idx]; }
                        if ( not defined $allMut_stat{$base}{lowCut}{$base} )  { $allMut_stat{$base}{lowCut}{$base} = $seq_mut{$seqID}{same}[$idx]; }
                        else  { $allMut_stat{$base}{lowCut}{$base} += $seq_mut{$seqID}{same}[$idx]; }
                    }
                    elsif ( $icSHAPE > $parameters{highCut} ) {
                        if ( not defined $seqMut_stat{$seqID}{$base}{highCut}{$base} )  { $seqMut_stat{$seqID}{$base}{highCut}{$base} = $seq_mut{$seqID}{same}[$idx]; }
                        else  { $seqMut_stat{$seqID}{$base}{highCut}{$base} += $seq_mut{$seqID}{same}[$idx]; }
                        if ( not defined $allMut_stat{$base}{highCut}{$base} )  { $allMut_stat{$base}{highCut}{$base} = $seq_mut{$seqID}{same}[$idx]; }
                        else  { $allMut_stat{$base}{highCut}{$base} += $seq_mut{$seqID}{same}[$idx]; }
                    }
                }
            }
            if ( defined $seq_mut{$seqID}{deletion}[$idx] ) {
                my $insD = () = $seq_mut{$seqID}{deletion}[$idx] =~ /;/gi; $insD++;

                if ( not defined $seqMut_stat{$seqID}{$base}{deletion} )  { $seqMut_stat{$seqID}{$base}{deletion} = $insD; }
                else { $seqMut_stat{$seqID}{$base}{deletion} += $insD; }
                if ( not defined $allMut_stat{$base}{deletion} )  { $allMut_stat{$base}{deletion} = $insD; }
                else { $allMut_stat{$base}{deletion} += $insD; }
                if ( $icSHAPE ne "NULL" ) {
                    if ( $icSHAPE < $parameters{lowCut} ) {
                        if ( not defined $seqMut_stat{$seqID}{$base}{lowCut}{deletion} )  { $seqMut_stat{$seqID}{$base}{lowCut}{deletion} = $insD; }
                        else { $seqMut_stat{$seqID}{$base}{lowCut}{deletion} += $insD; }
                        if ( not defined $allMut_stat{$base}{lowCut}{deletion} )  { $allMut_stat{$base}{lowCut}{deletion} = $insD; }
                        else { $allMut_stat{$base}{lowCut}{deletion} += $insD; }
                    }
                    elsif ( $icSHAPE > $parameters{highCut} ) {
                        if ( not defined $seqMut_stat{$seqID}{$base}{highCut}{deletion} )  { $seqMut_stat{$seqID}{$base}{highCut}{deletion} = $insD; }
                        else { $seqMut_stat{$seqID}{$base}{highCut}{deletion} += $insD; }
                        if ( not defined $allMut_stat{$base}{highCut}{deletion} )  { $allMut_stat{$base}{highCut}{deletion} = $insD; }
                        else { $allMut_stat{$base}{highCut}{deletion} += $insD; }
                    }
                }
            }
            if ( defined $seq_mut{$seqID}{mutation}[$idx] ) {
                my $insA = () = $seq_mut{$seqID}{mutation}[$idx] =~ /A/gi;
                my $insT = () = $seq_mut{$seqID}{mutation}[$idx] =~ /T/gi;
                my $insG = () = $seq_mut{$seqID}{mutation}[$idx] =~ /G/gi;
                my $insC = () = $seq_mut{$seqID}{mutation}[$idx] =~ /C/gi;

                if ( not defined $seqMut_stat{$seqID}{$base}{A} )  { $seqMut_stat{$seqID}{$base}{A} = $insA; }
                else { $seqMut_stat{$seqID}{$base}{A} += $insA; }
                if ( not defined $seqMut_stat{$seqID}{$base}{T} )  { $seqMut_stat{$seqID}{$base}{T} = $insT; }
                else { $seqMut_stat{$seqID}{$base}{T} += $insT; }
                if ( not defined $seqMut_stat{$seqID}{$base}{G} )  { $seqMut_stat{$seqID}{$base}{G} = $insG; }
                else { $seqMut_stat{$seqID}{$base}{G} += $insG; }
                if ( not defined $seqMut_stat{$seqID}{$base}{C} )  { $seqMut_stat{$seqID}{$base}{C} = $insC; }
                else { $seqMut_stat{$seqID}{$base}{C} += $insC; }

                if ( not defined $allMut_stat{$base}{A} )  { $allMut_stat{$base}{A} = $insA; }
                else { $allMut_stat{$base}{A} += $insA; }
                if ( not defined $allMut_stat{$base}{T} )  { $allMut_stat{$base}{T} = $insT; }
                else { $allMut_stat{$base}{T} += $insT; }
                if ( not defined $allMut_stat{$base}{G} )  { $allMut_stat{$base}{G} = $insG; }
                else { $allMut_stat{$base}{G} += $insG; }
                if ( not defined $allMut_stat{$base}{C} )  { $allMut_stat{$base}{C} = $insC; }
                else { $allMut_stat{$base}{C} += $insC; }

                if ( $icSHAPE ne "NULL" ) {
                    if ( $icSHAPE < $parameters{lowCut} ) {
                        if ( not defined $seqMut_stat{$seqID}{$base}{lowCut}{A} )  { $seqMut_stat{$seqID}{$base}{lowCut}{A} = $insA; }
                        else { $seqMut_stat{$seqID}{$base}{lowCut}{A} += $insA; }
                        if ( not defined $seqMut_stat{$seqID}{$base}{lowCut}{T} )  { $seqMut_stat{$seqID}{$base}{lowCut}{T} = $insT; }
                        else { $seqMut_stat{$seqID}{$base}{lowCut}{T} += $insT; }
                        if ( not defined $seqMut_stat{$seqID}{$base}{lowCut}{G} )  { $seqMut_stat{$seqID}{$base}{lowCut}{G} = $insG; }
                        else { $seqMut_stat{$seqID}{$base}{lowCut}{G} += $insG; }
                        if ( not defined $seqMut_stat{$seqID}{$base}{lowCut}{C} )  { $seqMut_stat{$seqID}{$base}{lowCut}{C} = $insC; }
                        else { $seqMut_stat{$seqID}{$base}{lowCut}{C} += $insC; }

                        if ( not defined $allMut_stat{$base}{lowCut}{A} )  { $allMut_stat{$base}{lowCut}{A} = $insA; }
                        else { $allMut_stat{$base}{lowCut}{A} += $insA; }
                        if ( not defined $allMut_stat{$base}{lowCut}{T} )  { $allMut_stat{$base}{lowCut}{T} = $insT; }
                        else { $allMut_stat{$base}{lowCut}{T} += $insT; }
                        if ( not defined $allMut_stat{$base}{lowCut}{G} )  { $allMut_stat{$base}{lowCut}{G} = $insG; }
                        else { $allMut_stat{$base}{lowCut}{G} += $insG; }
                        if ( not defined $allMut_stat{$base}{lowCut}{C} )  { $allMut_stat{$base}{lowCut}{C} = $insC; }
                        else { $allMut_stat{$base}{lowCut}{C} += $insC; }
                    }
                    elsif ( $icSHAPE > $parameters{highCut} ) {
                        if ( not defined $seqMut_stat{$seqID}{$base}{highCut}{A} )  { $seqMut_stat{$seqID}{$base}{highCut}{A} = $insA; }
                        else { $seqMut_stat{$seqID}{$base}{highCut}{A} += $insA; }
                        if ( not defined $seqMut_stat{$seqID}{$base}{highCut}{T} )  { $seqMut_stat{$seqID}{$base}{highCut}{T} = $insT; }
                        else { $seqMut_stat{$seqID}{$base}{highCut}{T} += $insT; }
                        if ( not defined $seqMut_stat{$seqID}{$base}{highCut}{G} )  { $seqMut_stat{$seqID}{$base}{highCut}{G} = $insG; }
                        else { $seqMut_stat{$seqID}{$base}{highCut}{G} += $insG; }
                        if ( not defined $seqMut_stat{$seqID}{$base}{highCut}{C} )  { $seqMut_stat{$seqID}{$base}{highCut}{C} = $insC; }
                        else { $seqMut_stat{$seqID}{$base}{highCut}{C} += $insC; }

                        if ( not defined $allMut_stat{$base}{highCut}{A} )  { $allMut_stat{$base}{highCut}{A} = $insA; }
                        else { $allMut_stat{$base}{highCut}{A} += $insA; }
                        if ( not defined $allMut_stat{$base}{highCut}{T} )  { $allMut_stat{$base}{highCut}{T} = $insT; }
                        else { $allMut_stat{$base}{highCut}{T} += $insT; }
                        if ( not defined $allMut_stat{$base}{highCut}{G} )  { $allMut_stat{$base}{highCut}{G} = $insG; }
                        else { $allMut_stat{$base}{highCut}{G} += $insG; }
                        if ( not defined $allMut_stat{$base}{highCut}{C} )  { $allMut_stat{$base}{highCut}{C} = $insC; }
                        else { $allMut_stat{$base}{highCut}{C} += $insC; }
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
            if ( defined $allMut_stat{$base}{$mut} ) {
                $total += $allMut_stat{$base}{$mut};
                print OUT "\t", $allMut_stat{$base}{$mut}; 
            }
            else { print OUT "\t0"; }
        }
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{$mut} ) { 
                if ( $base eq $mut ) { $percentage = "-" ; }
                else { $percentage = sprintf ( "%.5f", $allMut_stat{$base}{$mut} / $total ); }
                print OUT "\t", $percentage; 
            }
            else { print OUT "\t0"; }
        }

        $total = 0;
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{lowCut}{$mut} ) {
                $total += $allMut_stat{$base}{lowCut}{$mut};
                print OUT "\t", $allMut_stat{$base}{lowCut}{$mut}; 
            }
            else { print OUT "\t0"; }
        }
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{lowCut}{$mut} ) { 
                if ( $base eq $mut ) { $percentage = "-" ; }
                else { $percentage = sprintf ( "%.5f", $allMut_stat{$base}{lowCut}{$mut} / $total ); }
                print OUT "\t", $percentage; 
            }
            else { print OUT "\t0"; }
        }

        $total = 0;
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{highCut}{$mut} ) {
                $total += $allMut_stat{$base}{highCut}{$mut};
                print OUT "\t", $allMut_stat{$base}{highCut}{$mut}; 
            }
            else { print OUT "\t0"; }
        }
        foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
            if ( defined $allMut_stat{$base}{highCut}{$mut} ) { 
                if ( $base eq $mut ) { $percentage = "-" ; }
                else { $percentage = sprintf ( "%.5f", $allMut_stat{$base}{highCut}{$mut} / $total ); }
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
                if ( defined $seqMut_stat{$seqID}{$base}{$mut} ) { 
                    $total += $seqMut_stat{$seqID}{$base}{$mut};
                    print OUT "\t", $seqMut_stat{$seqID}{$base}{$mut}; 
                }
                else { print OUT "\t0"; }
            }
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $seqMut_stat{$seqID}{$base}{$mut} ) { 
                    if ( $base eq $mut ) { $percentage = "-" ; }
                    else { $percentage = sprintf ( "%.5f", $seqMut_stat{$seqID}{$base}{$mut} / $total ); }
                    print OUT "\t", $percentage; 
                }
                else { print OUT "\t0"; }
            }

            $total = 0;
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $seqMut_stat{$seqID}{$base}{lowCut}{$mut} ) { 
                    $total += $seqMut_stat{$seqID}{$base}{lowCut}{$mut};
                    print OUT "\t", $seqMut_stat{$seqID}{$base}{lowCut}{$mut}; 
                }
                else { print OUT "\t0"; }
            }
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $seqMut_stat{$seqID}{$base}{lowCut}{$mut} ) { 
                    if ( $base eq $mut ) { $percentage = "-" ; }
                    else { $percentage = sprintf ( "%.5f", $seqMut_stat{$seqID}{$base}{lowCut}{$mut} / $total ); }
                    print OUT "\t", $percentage; 
                }
                else { print OUT "\t0"; }
            }

            $total = 0;
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $seqMut_stat{$seqID}{$base}{highCut}{$mut} ) { 
                    $total += $seqMut_stat{$seqID}{$base}{highCut}{$mut};
                    print OUT "\t", $seqMut_stat{$seqID}{$base}{highCut}{$mut}; 
                }
                else { print OUT "\t0"; }
            }
            foreach my $mut ( "A", "T", "G", "C", "deletion" ) {
                if ( defined $seqMut_stat{$seqID}{$base}{highCut}{$mut} ) { 
                    if ( $base eq $mut ) { $percentage = "-" ; }
                    else { $percentage = sprintf ( "%.5f", $seqMut_stat{$seqID}{$base}{highCut}{$mut} / $total ); }
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
    my ( $seqID, $pos, $op ) = @_;

    my $match = 0;
    my $offset = 0;
    my $deletion = uc ( substr ( $op, 1 ) );
    my $substr = uc ( substr ( $seq{$seqID}, $pos-length($deletion)-1, 3*length($deletion) ) );
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
        for ( my $idx = 0; $idx < scalar ( @match ); $idx++ ) { if ( $match[$idx] eq "M" ) { if ( $matchSize[$idx] > $largestM ) { $largestM = $matchSize[$idx]; } } }

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

sub collectAlignInfoCIGAR
{
    my $seqID = shift;
    my $pos = shift;
    my $ref_match = shift;
    my $ref_matchSize = shift;
    my $readSeq = shift;
    my $cigar = shift;

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
            push ( @insertions, $refPos+1 );
            push ( @insertionSize, $ref_matchSize->[$idx] );
        }
        elsif ( $ref_match->[$idx] eq "P" ) { }
        elsif ( $ref_match->[$idx] eq "D" ) {
            $refPos += $ref_matchSize->[$idx];     
            # $seq_mut{$seqID}{"deletion"}[$refPos-1] = _append ( $seq_mut{$seqID}{"deletion"}[$refPos-1], $ref_matchSize->[$idx] ); # subject to realignment
        }
        elsif ( $ref_match->[$idx] eq "X" ) {
            $seq_mut{$seqID}{"mismatch"}[$refPos] = _append ( $seq_mut{$seqID}{"mismatch"}[$refPos], substr ( $readSeq, $readPos, $ref_matchSize->[$idx] ) ); ## not tested
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
    my $seqID = shift;
    my $pos = shift;
    my $headSoftClip = shift;
    my $ref_op = shift;
    my $ref_opSize = shift;
    my $readSeq = shift;
    my $ref_insertions = shift;
    my $ref_insertionSize = shift;

    my $refPos = $pos - 1;               ## point to 0-base index 
    my $readPos = $headSoftClip;         ## the same

    if ( $ref_op->[0] ) {
        if ( $ref_op->[0] =~ /^\^/ ) {
            ## removal of ambiguity
            my $delLen = length ( $ref_op->[0] ) - 1;
            for ( my $idxPos = 0; $idxPos < $delLen-1; $idxPos++ ) { $seq_mut{$seqID}{"same"}[$refPos+$idxPos]++; }
            $refPos += $delLen;
            if ( ( defined $ref_op->[1] ) and ( not defined $ref_opSize->[1] ) ) { print STDERR "error!\n"; return -1; }

            if ( ( defined $ref_op->[1] ) and ( $ref_op->[1] ) and ( not $ref_opSize->[0] ) ) { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            else {
                my $matchInLocal = _localReAlignment ( $seqID, $refPos-$delLen, $ref_op->[0] );
                print "multiple matching? - $matchInLocal\n" if ( $_debug );
                if ( $matchInLocal == 1 ) {
                    print "\tmatch only once, no realignment needed.\n" if ( $_debug );
                    $seq_mut{$seqID}{"deletion"}[$refPos-1] = _append ( $seq_mut{$seqID}{"deletion"}[$refPos-1], $ref_op->[0] );
                    ## you want to label the deletion even at the last base of deleted fragment
                }
                else { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            }
        }
        elsif ( $ref_op->[0] =~ /[A-Z]/ ) {
            $refPos++; $readPos++;
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            if ( ( defined $ref_op->[1] ) and ( not defined $ref_opSize->[1] ) ) { print STDERR "error!\n"; return -1; }

            if ( ( defined $ref_op->[1] ) and ( $ref_op->[1] ) and ( not $ref_opSize->[1] ) ) { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            else { $seq_mut{$seqID}{"mutation"}[$refPos-1] = _append ( $seq_mut{$seqID}{"mutation"}[$refPos-1], substr ( $readSeq, $readPos, 1 ) ); }
        }

        $ref_op->[0] = "";
        shift @{$ref_opSize};
    }

    for ( my $idx = 0; $idx <= scalar ( @{$ref_op} ); $idx++ ) {
        ## make sure the order is not affected by the real operation strings!!!!!!
        if ( ( not defined $ref_op->[$idx] ) or ( not $ref_op->[$idx] ) ) {
        }
        elsif ( $ref_op->[$idx] =~ /^\^/ ) {
            ## removal of ambiguity
            my $delLen = length ( $ref_op->[$idx] ) - 1;
            for ( my $idxPos = 0; $idxPos < $delLen-1; $idxPos++ ) { $seq_mut{$seqID}{"same"}[$refPos+$idxPos]++; }
            $refPos += $delLen;
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            if ( ( defined $ref_op->[$idx+1] ) and ( $ref_op->[$idx+1] ) and ( not $ref_opSize->[$idx] ) ) { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            else {
                my $matchInLocal = _localReAlignment ( $seqID, $refPos-$delLen, $ref_op->[$idx] );
                print "multiple matching? - $matchInLocal\n" if ( $_debug );
                if ( $matchInLocal == 1 ) {
                    print "\tmatch only once, no realignment needed.\n" if ( $_debug );
                    $seq_mut{$seqID}{"deletion"}[$refPos-1] = _append ( $seq_mut{$seqID}{"deletion"}[$refPos-1], $ref_op->[$idx] );
                    ## you want to label the deletion even at the last base of deleted fragment
                }
                else { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            }
        }
        elsif ( $ref_op->[$idx] =~ /[A-Z]/ ) {
            $refPos++; $readPos++;
            ## first check whether there is a following mutational event and whether it is immediately after this one. if true, skip this one
            if ( ( defined $ref_op->[$idx+1] ) and ( $ref_op->[$idx+1] ) and ( not $ref_opSize->[$idx] ) ) { $seq_mut{$seqID}{"same"}[$refPos-1]++; }
            else { $seq_mut{$seqID}{"mutation"}[$refPos-1] = _append ( $seq_mut{$seqID}{"mutation"}[$refPos-1], substr ( $readSeq, $readPos, 1 ) ); }
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
