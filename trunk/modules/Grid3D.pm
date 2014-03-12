# Copyright 2012 Thomas H. Schmidt
#
# This file is part of squaredance.
#
# squaredance is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# squaredance is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with squaredance; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

package Grid3D;

use strict;
use warnings;
use pte;
use FileIO::gro;

#require Exporter;
our $VERSION = 1.0;
#our @ISA     = qw(Exporter);
#our @EXPORT  = qw(atoms2Grid getGridRef grid2GroFile setGridDelta);


my @gridOcc;
my @gridInt;
my @vdwSurf;

my @gridXMin;
my @gridXMax;
my @gridYMin;
my @gridYMax;
my $gridZMin = 1000;
my $gridZMax = 0;

my $gridDeltaX = 0.5;
my $gridDeltaY = 0.5;
my $gridDeltaZ = 0.5;

my $coordShift = 2; # Shift all coordinates to avoid negative grid points (PBC). FIND A BETTER SOLUTION FOR THIS!



sub setGridDelta {
    return 0 unless defined $_[0];
    $gridDeltaX = $_[0];
    $gridDeltaY = defined $_[1] ? $_[1] : $_[0];
    $gridDeltaZ = defined $_[2] ? $_[2] : $_[1];
    return 1;
}

sub getGridOccRef { return \@gridOcc; }

sub getGridIntRef { return \@gridInt; }

sub getGridVdwSurfRef { return \@vdwSurf; }



sub getMin {
    my $min = $_[$gridZMin];
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        next unless defined $_[$z];
        $min = $_[$z] if $_[$z] < $min;
    }
    return $min;
}



sub getMax {
    my $max = $_[$gridZMax];
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        next unless defined $_[$z];
        $max = $_[$z] if $_[$z] > $max;
    }
    return $max;
}



sub atoms2Grid {
    my $atomsIdsRef  = shift;
    my $coordDataRef = shift;
    my $gridOccVal   = shift;

    my $atomId   = 0;

    print "  ---------------------------------\n  Mapping atoms to the grid\r";

    foreach (@{$atomsIdsRef}) {
        next unless $$coordDataRef[$_]{'cooZ'};

        my $element = substr($$coordDataRef[$_]{'atomName'}, 0, 1);
#        my $element = $$coordDataRef[$_]{'atomName'};
        my $radius  = PTE::getRadius($element)*1.4;
        my $radius2 = $radius * $radius;


        my $tmpGridX = sprintf("%d", round(($$coordDataRef[$_]{'cooX'} + $coordShift) / $gridDeltaX, 1));
        my $tmpGridY = sprintf("%d", round(($$coordDataRef[$_]{'cooY'} + $coordShift) / $gridDeltaY, 1));
        my $tmpGridZ = sprintf("%d", round(($$coordDataRef[$_]{'cooZ'} + $coordShift) / $gridDeltaZ, 1));
#        my $subrangeX = $radius > $gridDeltaX ? sprintf("%d", round($radius / $gridDeltaX, 1)) : 0;
        my $subrangeX = sprintf("%d", round($radius / $gridDeltaX, 1));
#        my $subrangeY = $radius > $gridDeltaY ? sprintf("%d", round($radius / $gridDeltaY, 1)) : 0;
        my $subrangeY = sprintf("%d", round($radius / $gridDeltaY, 1));
#        my $subrangeZ = $radius > $gridDeltaZ ? sprintf("%d", round($radius / $gridDeltaZ, 1)) : 0;
        my $subrangeZ = sprintf("%d", round($radius / $gridDeltaZ, 1));

        for (my $z=($tmpGridZ-$subrangeZ); $z<=($tmpGridZ+$subrangeZ); $z++) {
            for (my $x=($tmpGridX-$subrangeX); $x<=($tmpGridX+$subrangeX); $x++) {
                for (my $y=($tmpGridY-$subrangeY); $y<=($tmpGridY+$subrangeY); $y++) {
                    my $dx = $$coordDataRef[$_]{'cooX'} + $coordShift - $x * $gridDeltaX;
                    my $dy = $$coordDataRef[$_]{'cooY'} + $coordShift - $y * $gridDeltaY;
                    next if ($dx*$dx + $dy*$dy) > $radius2;
                    my $dz = $$coordDataRef[$_]{'cooZ'} + $coordShift - $z * $gridDeltaZ;
                    next if ($dx*$dx + $dz*$dz) > $radius2;
                    next if ($dy*$dy + $dz*$dz) > $radius2;

                    $gridOcc[$z][$x][$y] = $gridOccVal ? $gridOccVal : $$coordDataRef[$_]{'resName'};

                    $gridXMin[$z] = $x if !defined $gridXMin[$z] || $x < $gridXMin[$z];
                    $gridXMax[$z] = $x if !defined $gridXMax[$z] || $x > $gridXMax[$z];
                    $gridYMin[$z] = $y if !defined $gridYMin[$z] || $y < $gridYMin[$z];
                    $gridYMax[$z] = $y if !defined $gridYMax[$z] || $y > $gridYMax[$z];
                    $gridZMin     = $z if $z < $gridZMin;
                    $gridZMax     = $z if $z > $gridZMax;
                }
            }
        }
        printf("  Mapping atoms to the grid: %d%%\r", ++$atomId*100/@{$atomsIdsRef}) if $main::verbose;
    }

    printf("  Mapping atoms to the grid: Finished\n  ---------------------------------\n\n") if $main::verbose;

    return 1;
}



sub getVdwSurf {
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        next unless $gridInt[$z];
        for (my $x=$gridXMin[$z]; $x<=$gridXMax[$z]; $x++) {
            for (my $y=$gridYMin[$z]; $y<=$gridYMax[$z]; $y++) {
                next unless $gridInt[$z][$x][$y];
                next if $gridInt[$z][$x][$y] eq 'INT';
                $vdwSurf[$z][$x+1][$y]   = $gridOcc[$z][$x+1][$y]   if ($gridOcc[$z][$x+1][$y]);
                $vdwSurf[$z][$x+1][$y+1] = $gridOcc[$z][$x+1][$y+1] if ($gridOcc[$z][$x+1][$y+1]);
                $vdwSurf[$z][$x][$y+1]   = $gridOcc[$z][$x][$y+1]   if ($gridOcc[$z][$x][$y+1]);
                $vdwSurf[$z][$x-1][$y+1] = $gridOcc[$z][$x-1][$y+1] if ($gridOcc[$z][$x-1][$y+1]);
                $vdwSurf[$z][$x-1][$y]   = $gridOcc[$z][$x-1][$y]   if ($gridOcc[$z][$x-1][$y]);
                $vdwSurf[$z][$x-1][$y-1] = $gridOcc[$z][$x-1][$y-1] if ($gridOcc[$z][$x-1][$y-1]);
                $vdwSurf[$z][$x][$y-1]   = $gridOcc[$z][$x][$y-1]   if ($gridOcc[$z][$x][$y-1]);
                $vdwSurf[$z][$x+1][$y-1] = $gridOcc[$z][$x+1][$y-1] if ($gridOcc[$z][$x+1][$y-1]);
            }
        }
    }
}


sub countVoxelsWithinZRange {
    my $gridRef     = shift;
    my $zLimitLower = shift;
    my $zLimitUpper = shift;
    my $valueRegex  = shift;

    my $zLimitLowerGrid = sprintf("%d", round(($zLimitLower + $coordShift) / $gridDeltaZ, 1));
    my $zLimitUpperGrid = sprintf("%d", round(($zLimitUpper + $coordShift) / $gridDeltaZ, 1));
    my $nVoxels = 0;
    my $nVoxelsTotal = 0;

    for (my $z=$zLimitLowerGrid; $z<=$zLimitUpperGrid; $z++) {
        next unless $$gridRef[$z];
        for (my $x=$gridXMin[$z]; $x<=$gridXMax[$z]; $x++) {
            next unless $$gridRef[$z][$x];
            for (my $y=$gridYMin[$z]; $y<=$gridYMax[$z]; $y++) {
                next unless $$gridRef[$z][$x][$y];
                if ($valueRegex) {
                    $nVoxels++ if $$gridRef[$z][$x][$y] =~ /$valueRegex/;
                }
                else {
                    $nVoxels++
                }
                $nVoxelsTotal++;
            }
        }
    }
    return $nVoxels*100/$nVoxelsTotal;
}



sub getSas {
    my $radius = shift;
    $radius    = 0.14 unless defined $radius; # This could also be 0 to washing out cavities with even small gaps to the exterior.
    my @xyProbeTemplate = createXYProbeTemplate($radius);

    my $nCavityAreas = 0;
    my $gridDelta2   = $gridDeltaX * $gridDeltaY;

    ### Plus imaginary space for surface voxels and cavity analysis ############
    extendGridXY();
    ############################################################################

    print "  ---------------------------------\n  Detecting internal cavities\r" if $main::verbose;

    ### Invert existing grid (define void voxels as cavity) ####################
    $nCavityAreas = invertGrid(\@gridOcc, \@gridInt, 'INT');
    ############################################################################


    ### Washing out the "external" cavities ####################################
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        for (my $i=0; $i<2; $i++) {
            my $foundExcl1 = 1;
            my $foundExcl2 = 1;
            for (my $x=$gridXMin[$z]; $x<=$gridXMax[$z]; $x++) { # W -> E.
                $nCavityAreas -= $foundExcl1 = washingOutYSas($gridOcc[$z], $gridInt[$z], $x, $z, $x-1, \@xyProbeTemplate);
            }

            for (my $x=$gridXMax[$z]; $x>=$gridXMin[$z]; $x--) { # E -> W.
                $nCavityAreas -= $foundExcl2 = washingOutYSas($gridOcc[$z], $gridInt[$z], $x, $z, $x+1, \@xyProbeTemplate);
            }
            $i = 0 if $foundExcl1 || $foundExcl2;
            printf("  Detecting internal cavities: %.4f nm^2 possible      \r", $nCavityAreas*$gridDelta2) if $main::verbose;
        }
    }
    ############################################################################

    printf("  Detecting internal cavities: Finished\n  ---------------------------------\n\n", $nCavityAreas*$gridDelta2) if $main::verbose;

    return ($nCavityAreas*$gridDelta2);
}



sub detectCavity {
    my $radius = shift;
    $radius    = 0.14 unless defined $radius; # This could also be 0 to washing out cavities with even small gaps to the exterior.
    my @xyProbeTemplate = createXYProbeTemplate($radius);

    my $nCavityAreas = 0;
    my $gridDelta2   = $gridDeltaX * $gridDeltaY;

    ### Plus imaginary space for surface voxels and cavity analysis ############
    extendGridXY(8);
    ############################################################################

    print "  ---------------------------------\n  Detecting internal cavities\r" if $main::verbose;

    ### Invert existing grid (define void voxels as cavity) ####################
    $nCavityAreas = invertGrid(\@gridOcc, \@gridInt, 'INT');
    ############################################################################


    ### Washing out the "external" cavities ####################################
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        for (my $i=0; $i<2; $i++) {
            printf("  Detecting internal cavities: %.4f nm^2 possible      \r", $nCavityAreas*$gridDelta2) if $main::verbose;
            my $foundExcl1 = 1;
            my $foundExcl2 = 1;
            for (my $x=$gridXMin[$z]; $x<=$gridXMax[$z]; $x++) { # W -> E.
                $nCavityAreas -= $foundExcl1 = washingOutY2($gridOcc[$z], $gridInt[$z], $x, $z, $x-1, \@xyProbeTemplate);
            }

            for (my $x=$gridXMax[$z]; $x>=$gridXMin[$z]; $x--) { # E -> W.
                $nCavityAreas -= $foundExcl2 = washingOutY2($gridOcc[$z], $gridInt[$z], $x, $z, $x+1, \@xyProbeTemplate);
            }
            $i = 0 if $foundExcl1 || $foundExcl2;
        }
    }
    ############################################################################

    printf("  Detecting internal cavities: Finished\n  ---------------------------------\n\n", $nCavityAreas*$gridDelta2) if $main::verbose;

#    ### Updating INT grid (separation of SUR voxels) ############################
#    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
#        
#    }
#    ############################################################################
    

    return ($nCavityAreas*$gridDelta2);
}



sub invertGrid {
    my $gridOccRef = shift;
    my $gridInvRef = shift;
    my $voxelName  = shift;
    my $nInvVoxels  = 0;
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        for (my $x=$gridXMin[$z]; $x<=$gridXMax[$z]; $x++) {
            for (my $y=$gridYMin[$z]; $y<=$gridYMax[$z]; $y++) {
                next if defined $$gridOccRef[$z][$x][$y];
                $$gridInvRef[$z][$x][$y] = $voxelName;
                $nInvVoxels++;
            }
        }
    }
    return $nInvVoxels;
}



sub extendGridXY {
    my $extLen = shift;
    $extLen    = 1 unless $extLen;
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        $gridXMin[$z]-=$extLen;
        $gridYMin[$z]-=$extLen;
        $gridXMax[$z]+=$extLen;
        $gridYMax[$z]+=$extLen;
    }
}



sub round {
    my ( $num, $prec ) = @_;
    return int( $num / $prec + 0.5 - ( $num < 0 ) ) * $prec;
}


my @probe;
my @xyProbeTemplate;
sub getProbeRef {
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        for (my $i=0; $i<@xyProbeTemplate; $i++) {
            $probe[$z][$gridXMin[$z]+$xyProbeTemplate[$i]{'x'}][$gridYMin[$z]+$xyProbeTemplate[$i]{'y'}] = 1;
        }
    }
    return \@probe;
}


sub createXYProbeTemplate {
    my $radius  = shift;
    return 0 unless $radius;
    my $radius2 = $radius * $radius;

    my $subrangeX = $radius > $gridDeltaX ? sprintf("%d", round($radius / $gridDeltaX, 1)) : 0;
    my $subrangeY = $radius > $gridDeltaY ? sprintf("%d", round($radius / $gridDeltaY, 1)) : 0;

    for (my $x=(0-$subrangeX); $x<=(0+$subrangeX); $x++) {
        for (my $y=(0-$subrangeY); $y<=(0+$subrangeY); $y++) {
            my $dx = $x * $gridDeltaX;
            my $dy = $y * $gridDeltaY;
            next if ($dx*$dx + $dy*$dy) > $radius2;
            my %tmp = ('x' => $x, 'y' => $y);
            push(@xyProbeTemplate, \%tmp);
        }
    }
    return @xyProbeTemplate;
}



sub washingOutY {
    my $gridOccSliceRef = shift;
    my $gridIntSliceRef = shift;
    my $x               = shift;
    my $z               = shift;
    my $neighbCellX     = shift;
    my $xyProbeTemplateRef = shift;

    my $delAreas        = 0;


    for (my $y=$gridYMin[$z]; $y<=$gridYMax[$z]; $y++) { # S -> N.
        next unless $$gridIntSliceRef[$x][$y]; # Next if this grid point is not a possible cavity.
        my $neighbCellY = $y - 1;
        if ((!$$gridIntSliceRef[$x][$neighbCellY] && !$$gridOccSliceRef[$x][$neighbCellY]) ||
            (!$$gridIntSliceRef[$neighbCellX][$y] && !$$gridOccSliceRef[$neighbCellX][$y])) { # If the neighbored cell is not protein or not a possible cavity...
            unless (searchOverlaps($gridOccSliceRef, $gridIntSliceRef, $x, $y, 5, 5, $z)) {
                $$gridIntSliceRef[$x][$y] = 0; # ...exclude this from a possible cavity.
                $delAreas++;
                next;
            }
        }
    }

    for (my $y=$gridYMax[$z]; $y>=$gridYMin[$z]; $y--) { # N -> S.
        next unless $$gridIntSliceRef[$x][$y];
        my $neighbCellY = $y + 1;
        if ((!$$gridIntSliceRef[$x][$neighbCellY] && !$$gridOccSliceRef[$x][$neighbCellY]) ||
            (!$$gridIntSliceRef[$neighbCellX][$y] && !$$gridOccSliceRef[$neighbCellX][$y])) {
            unless (searchOverlaps($gridOccSliceRef, $gridIntSliceRef, $x, $y, 5, 5, $z)) {
                $$gridIntSliceRef[$x][$y] = 0;
                $delAreas++;
                next;
            }
        }
    }

    return $delAreas;
}



sub washingOutY2 {
    my $gridOccSliceRef = shift;
    my $gridIntSliceRef = shift;
    my $x               = shift;
    my $z               = shift;
    my $neighbCellX     = shift;
    my $xyProbeTemplateRef = shift;

    my $delAreas        = 0;

    for (my $y=$gridYMin[$z]; $y<=$gridYMax[$z]; $y++) { # S -> N.
        next unless $$gridIntSliceRef[$x][$y]; # Next if this grid point is not a possible cavity.
        my $neighbCellY = $y - 1;
        if ((!$$gridIntSliceRef[$x][$neighbCellY] && !$$gridOccSliceRef[$x][$neighbCellY]) ||
            (!$$gridIntSliceRef[$neighbCellX][$y] && !$$gridOccSliceRef[$neighbCellX][$y])) { # If the neighbored cell is not protein or not a possible cavity...
            unless (checkXyProbeOverlaps($gridOccSliceRef, $xyProbeTemplateRef, $x, $y)) {
                $delAreas += delXyProbeOverlaps($gridIntSliceRef, $xyProbeTemplateRef, $x, $y);
                $$gridIntSliceRef[$x][$y] = 0;
                next;
            }
        }
    }

    for (my $y=$gridYMax[$z]; $y>=$gridYMin[$z]; $y--) { # N -> S.
        next unless $$gridIntSliceRef[$x][$y];
        my $neighbCellY = $y + 1;
        if ((!$$gridIntSliceRef[$x][$neighbCellY] && !$$gridOccSliceRef[$x][$neighbCellY]) ||
            (!$$gridIntSliceRef[$neighbCellX][$y] && !$$gridOccSliceRef[$neighbCellX][$y])) {
            unless (checkXyProbeOverlaps($gridOccSliceRef, $xyProbeTemplateRef, $x, $y)) {
                $delAreas += delXyProbeOverlaps($gridIntSliceRef, $xyProbeTemplateRef, $x, $y);
                $$gridIntSliceRef[$x][$y] = 0;
                next;
            }
        }
    }

    return $delAreas;
}



sub washingOutYSas {
    my $gridOccSliceRef = shift;
    my $gridIntSliceRef = shift;
    my $x               = shift;
    my $z               = shift;
    my $neighbCellX     = shift;
    my $xyProbeTemplateRef = shift;

    my $delAreas        = 0;

    for (my $y=$gridYMin[$z]; $y<=$gridYMax[$z]; $y++) { # S -> N.
        next unless $$gridIntSliceRef[$x][$y]; # Next if this grid point is not a possible cavity.
        my $neighbCellY = $y - 1;
#        print "\nWill wash out this grid point, if its neighbors are not defined " . $$gridOccSliceRef[$x][$neighbCellY] . ", " . $$gridIntSliceRef[$x][$neighbCellY] . "\n" if $x == 75 && $y == 120 && $z == 213;
#        unless (($$gridIntSliceRef[$x][$neighbCellY] && $$gridOccSliceRef[$x][$neighbCellY]) || $$gridIntSliceRef[$x][$neighbCellY] eq 'INT') { # If the neighbored cell is not protein or not a possible cavity...
        if ((!$$gridIntSliceRef[$x][$neighbCellY] && !$$gridOccSliceRef[$x][$neighbCellY]) ||
            (!$$gridIntSliceRef[$neighbCellX][$y] && !$$gridOccSliceRef[$neighbCellX][$y])) { # If the neighbored cell is not protein or not a possible cavity...
            unless (checkXyProbeOverlaps($gridOccSliceRef, $xyProbeTemplateRef, $x, $y)) {
            #searchOverlaps($gridOccSliceRef, $gridIntSliceRef, $x, $y, $maxRadiusGridX, $maxRadiusGridY)) {
                $$gridIntSliceRef[$x][$y] = 0; # ...exclude this from a possible cavity.
                $delAreas++;
                next;
            }
        }
#        if (!$$gridIntSliceRef[$neighbCellX][$y] && !$$gridOccSliceRef[$neighbCellX][$y]) {
#            unless (checkXyProbeOverlaps($gridOccSliceRef, $xyProbeTemplateRef, $x, $y)) {
#            #unless (searchOverlaps($gridOccSliceRef, $gridIntSliceRef, $x, $y, $maxRadiusGridX, $maxRadiusGridY)) {
#                $$gridIntSliceRef[$x][$y] = 0;
#                $delAreas++;
#                next;
#            }
#        }
    }

    for (my $y=$gridYMax[$z]; $y>=$gridYMin[$z]; $y--) { # N -> S.
        next unless $$gridIntSliceRef[$x][$y];
        my $neighbCellY = $y + 1;
        if ((!$$gridIntSliceRef[$x][$neighbCellY] && !$$gridOccSliceRef[$x][$neighbCellY]) ||
            (!$$gridIntSliceRef[$neighbCellX][$y] && !$$gridOccSliceRef[$neighbCellX][$y])) {
            unless (checkXyProbeOverlaps($gridOccSliceRef, $xyProbeTemplateRef, $x, $y)) {
            #unless (searchOverlaps($gridOccSliceRef, $gridIntSliceRef, $x, $y, $maxRadiusGridX, $maxRadiusGridY)) {
                $$gridIntSliceRef[$x][$y] = 0;
                $delAreas++;
                next;
            }
        }
#        if (!$$gridIntSliceRef[$neighbCellX][$y] && !$$gridOccSliceRef[$neighbCellX][$y]) {
#            unless (checkXyProbeOverlaps($gridOccSliceRef, $xyProbeTemplateRef, $x, $y)) {
#            #unless (searchOverlaps($gridOccSliceRef, $gridIntSliceRef, $x, $y, $maxRadiusGridX, $maxRadiusGridY)) {
#                $$gridIntSliceRef[$x][$y] = 0;
#                $delAreas++;
#                next;
#            }
#        }
    }

    return $delAreas;
}



sub checkXyProbeOverlaps {
    my $gridOccSliceRef    = shift;
    my $xyProbeTemplateRef = shift;
    my $x = shift;
    my $y = shift;
    my $nReqOverlaps = shift || 1;
    my $nOverlaps = 0;
    foreach (@$xyProbeTemplateRef) {
        if ($$gridOccSliceRef[$x+$$_{'x'}][$y+$$_{'y'}]) {
            $nOverlaps++;
            return 1 if $nOverlaps == $nReqOverlaps;
        }
    }
    return 0;
}



sub delXyProbeOverlaps {
    my $gridIntSliceRef    = shift;
    my $xyProbeTemplateRef = shift;
    my $x = shift;
    my $y = shift;
    my $nVoxels = 0;
    foreach (@$xyProbeTemplateRef) {
        next unless $$gridIntSliceRef[$x+$$_{'x'}][$y+$$_{'y'}];
        $$gridIntSliceRef[$x+$$_{'x'}][$y+$$_{'y'}] = 'SAS';
        $nVoxels++;
    }
    return $nVoxels;
}



sub searchOverlaps {
    my $gridOccSliceRef = shift;
    my $gridIntSliceRef = shift;
    my $x               = shift;
    my $y               = shift;
    my $maxRadiusGridX  = shift;
    my $maxRadiusGridY  = shift;
    my $z               = shift;

    my $xMin = $x-$maxRadiusGridX;
    my $xMax = $x+$maxRadiusGridX;
    my $yMin = $y-$maxRadiusGridY;
    my $yMax = $y+$maxRadiusGridY;


    ### Vertical ###############################################################
    my $protLimitMin = 0;
    my $protLimitMax = 0;
    my $nCellsMin    = 0;
    my $nCellsMax    = 0;
    my $tmpY         = $y;
    while (!$protLimitMin && $tmpY >= $yMin) {
        $protLimitMin = 1 if $$gridOccSliceRef[$x][$tmpY--];
        $nCellsMin--;
    }
    $tmpY = $y;
    while (!$protLimitMax && $tmpY <= $yMax) {
        $protLimitMax = 1 if $$gridOccSliceRef[$x][$tmpY++];
        $nCellsMax++;
    }
    if ($protLimitMin && $protLimitMax) {
        printf("  Cell (%f,%f) is limited by protein, forming a cavity with cells (%f,%f-%f)\n", ($x*$gridDeltaX-$coordShift)*10, ($y*$gridDeltaY-$coordShift)*10, ($x*$gridDeltaX-$coordShift)*10, (($y-$nCellsMin)*$gridDeltaY-$coordShift)*10, (($y+$nCellsMax)*$gridDeltaY-$coordShift)*10, ($y*$gridDeltaY-$coordShift)*10) if $z == 53;
        for (my $cavY=$nCellsMin; $cavY<=$nCellsMax; $cavY++) {
            $$gridIntSliceRef[$x][$cavY] = 'CAV';
        }
        return 1;
    }
    ############################################################################
    return 0;


    ### Horizontal #############################################################
    $protLimitMin = 0;
    $protLimitMax = 0;
    $nCellsMin    = 0;
    $nCellsMax    = 0;
    my $tmpX      = $x;
    while (!$protLimitMin && $tmpX >= $xMin) {
        $protLimitMin = 1 if $$gridOccSliceRef[$tmpX--][$y];
        $nCellsMin--;
    }
    $tmpX = $x;
    while (!$protLimitMax && $tmpX <= $xMax) {
        $protLimitMax = 1 if $$gridOccSliceRef[$tmpX++][$y];
        $nCellsMax++;
    }
    if ($protLimitMin && $protLimitMax) {
        for (my $cavX=$nCellsMin; $cavX<=$nCellsMax; $cavX++) {
            $$gridIntSliceRef[$cavX][$y] = 'CAV';
        }
        return 1;
    }
    ############################################################################

    return 0;
#
#
#    if ($$gridOccSliceRef[$x][$yMin]) {
#        my $tmpY = $y;
#        while ($tmpY < $yMax) {
#            $tmpY = $yMax if !$$gridIntSliceRef[$x][$tmpY] && !$$gridOccSliceRef[$x][$tmpY];
#            return 1 if $$gridOccSliceRef[$x][$tmpY++];
#        }
#    }
#    if ($$gridOccSliceRef[$x][$yMax]) {
#        my $tmpY = $y;
#        while ($tmpY > $yMin) {
#            $tmpY = $yMin if !$$gridIntSliceRef[$x][$tmpY] && !$$gridOccSliceRef[$x][$tmpY];
#            return 1 if $$gridOccSliceRef[$x][$tmpY--];
#        }
#    }
#    if ($$gridOccSliceRef[$xMin][$y]) {
#        my $tmpX = $x;
#        while ($tmpX < $xMax) {
#            $tmpX = $xMax if !$$gridIntSliceRef[$tmpX][$y] && !$$gridOccSliceRef[$tmpX][$y];
#            return 1 if $$gridOccSliceRef[$tmpX++][$y];
#        }
#    }
#    if ($$gridOccSliceRef[$xMax][$y]) {
#        my $tmpX = $x;
#        while ($tmpX > $xMin) {
#            $tmpX = $xMin if !$$gridIntSliceRef[$tmpX][$y] && !$$gridOccSliceRef[$tmpX][$y];
#            return 1 if $$gridOccSliceRef[$tmpX--][$y];
#        }
#    }
#    if ($$gridOccSliceRef[$xMin][$yMin]) {
#        my $tmpX = $x;
#        my $tmpY = $y;
#        while ($tmpX < $xMax && $tmpY < $yMax) {
#            return 1 if $$gridOccSliceRef[$tmpX++][$tmpY++];
#        }
#    }
#    if ($$gridOccSliceRef[$xMin][$yMax]) {
#        my $tmpX = $x;
#        my $tmpY = $y;
#        while ($tmpX < $xMax && $tmpY > $yMin) {
#            return 1 if $$gridOccSliceRef[$tmpX++][$tmpY--];
#        }
#    }
#    return 0;
}



sub grid2GroFile {
    my $gridRef      = shift;
    my $gridOutFile  = shift;

    my $voxId        = 0;
    my %gridGroData;

    my $gridXMinAbs = getMin(@gridXMin);
    my $gridXMaxAbs = getMax(@gridXMax);
    my $gridYMinAbs = getMin(@gridYMin);
    my $gridYMaxAbs = getMax(@gridYMax);

    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        for (my $x=$gridXMinAbs; $x<=$gridXMaxAbs; $x++) {
            for (my $y=$gridYMinAbs; $y<=$gridYMaxAbs; $y++) {
                next unless $$gridRef[$z][$x][$y];
                my $resName = $$gridRef[$z][$x][$y] ? $$gridRef[$z][$x][$y] : 'MOX';
                my %tmpCoords = ('cooX' => $x * $gridDeltaX - $coordShift,
                                 'cooY' => $y * $gridDeltaY - $coordShift,
                                 'cooZ' => $z * $gridDeltaZ - $coordShift);
                push(@{$gridGroData{'atoms'}}, setVoxel(++$voxId, $resName, $$gridRef[$z][$x][$y], $voxId, \%tmpCoords));
            }
        }
    }
    $gridGroData{'title'}  = sprintf("Grid stepwidth x = %f, y = %f, z = %f", $gridDeltaX, $gridDeltaY, $gridDeltaZ);
    $gridGroData{'box'}{'cooX'} = $gridXMaxAbs * $gridDeltaX;
    $gridGroData{'box'}{'cooY'} = $gridYMaxAbs * $gridDeltaY;
    $gridGroData{'box'}{'cooZ'} = $gridZMax * $gridDeltaZ;
    $gridGroData{'atoms'}  = GRO::renumAtoms(\@{$gridGroData{'atoms'}});
    $gridGroData{'nAtoms'} = scalar(@{$gridGroData{'atoms'}}) - 1;
    GRO::writeGro($gridOutFile, \%gridGroData);
}



sub setVoxel {
    my $residue   = shift;
    my $resName   = shift;
    my $atomName  = shift;
    my $serial    = shift;
    my $coordsRef = shift;
    my %voxel = ('residue'  => $residue,
                 'resId'    => $residue,
                 'resName'  => $resName,
                 'atomName' => $atomName,
                 'serial'   => $serial,
                 'cooX'     => $$coordsRef{'cooX'},
                 'cooY'     => $$coordsRef{'cooY'},
                 'cooZ'     => $$coordsRef{'cooZ'});
    return \%voxel;
}



sub writeXyProfile {
    my $xyProfile = shift;
    my $gridRef   = shift;
    my $gridsizeX = shift;
    my $gridsizeY = shift;
    my $gridDelta = shift;

    printf("      Writing out protein xy grid profile: 0%%\r") if $main::verbose;
    open(XYPROFILE, ">$xyProfile") || die "ERROR: Cannot open profile output file \"$xyProfile\": $!\n";
    for (my $x=0; $x<=$gridsizeX; $x++) {
        my $tmpX = $x*$gridDelta; # Performance.
        for (my $y=0; $y<=$gridsizeY; $y++) {
            printf(XYPROFILE "%f %f %d\n", $tmpX, $y*$gridDelta, $$gridRef[$x][$y]) if $$gridRef[$x][$y];
        }
        printf("      Writing out protein xy grid profile: %d%%\r", 100*$x/$gridsizeX) if $main::verbose;
    }
    close XYPROFILE;
    printf("      Writing out protein xy grid profile: 100%%\n") if $main::verbose;
}



#sub connectSubunits {
#    my $gridRef       = shift;
#    my $gridsizeX     = shift;
#    my $gridsizeY     = shift;
#    my $gridDelta     = shift;
#    my $ndxDataRef    = shift;
#    my $suGroupIdsRef = shift;
#    my $coordsRef     = shift;
#    my $zMin          = shift;
#    my $zMax          = shift;
#
#    my @gridSuGeoCenter;
#
#
#    print "      Detecting subunit connections: 0.0000 nm^2 possible\r" if $main::verbose;
#
#    ### Calculate center of mass of each protein subunit in that slice #########
#    foreach my $suGroupId (@{$suGroupIdsRef}) {
#        my %tmpSum = ('cooX' => 0, 'cooY' => 0);
#        my $nAtoms = 0;
#        foreach (@{$$ndxDataRef[$suGroupId]{'atoms'}}) {
#            next unless $$coordsRef[$_]{'cooZ'};
#            next if $$coordsRef[$_]{'cooZ'} > $zMax;
#            next if $$coordsRef[$_]{'cooZ'} < $zMin;
#            $tmpSum{'cooX'} += $$coordsRef[$_]{'cooX'};
#            $tmpSum{'cooY'} += $$coordsRef[$_]{'cooY'};
#            $nAtoms++;
#        }
#        next unless $nAtoms;
#        $tmpSum{'cooX'} = int(($tmpSum{'cooX'}/$nAtoms)/$gridDelta);
#        $tmpSum{'cooY'} = int(($tmpSum{'cooY'}/$nAtoms)/$gridDelta);
#        push(@gridSuGeoCenter, \%tmpSum);
#    }
#    ############################################################################
#
#
#    ### Print the geometrical center of each subunit ###########################
#    for (my $i=0; $i<@gridSuGeoCenter; $i++) {
#        getGridSlope($gridRef, $gridSuGeoCenter[$i-1], $gridSuGeoCenter[$i]);
#    }
#    ############################################################################
#}



#sub getGridSlope {
#    my $gridRef   = shift;
#    my $vecARef   = shift;
#    my $vecBRef   = shift;
#
#    my $slope = ($$vecBRef{'cooY'} - $$vecARef{'cooY'}) / ($$vecBRef{'cooX'} - $$vecARef{'cooX'});
#    my $xMin = $$vecARef{'cooX'} < $$vecBRef{'cooX'} ? $$vecARef{'cooX'} : $$vecBRef{'cooX'};
#    my $xMax = $$vecARef{'cooX'} > $$vecBRef{'cooX'} ? $$vecARef{'cooX'} : $$vecBRef{'cooX'};
#
#    my $y = $$vecARef{'cooY'};
#    for (my $x=$xMin; $x<=$xMax; $x++) {
#        $y += $slope;
#        my $intY = int($y);
#        next if $$gridRef[$x][$intY];
#        $$gridRef[$x][$intY] = 3;
#
#        $$gridRef[$x+1][$intY] = 3;
#        $$gridRef[$x-1][$intY] = 3;
#        $$gridRef[$x][$intY+1] = 3;
#        $$gridRef[$x][$intY-1] = 3;
#        $$gridRef[$x+1][$intY+1] = 3;
#        $$gridRef[$x+1][$intY-1] = 3;
#        $$gridRef[$x-1][$intY+1] = 3;
#        $$gridRef[$x-1][$intY-1] = 3;
#    }
#}



#sub ini3DGrid {
#    my $boxRef     = shift;
#    my $gridDeltaZ = shift;
#    my $gridDeltaX = shift;
#    my $gridDeltaY = shift;
#
#    $gridDeltaZ = 0.5 unless $gridDeltaZ;
#    $gridDeltaX = $gridDeltaZ unless $gridDeltaX;
#    $gridDeltaY = $gridDeltaX unless $gridDeltaY;
#
#    my $gridsizeX = sprintf("%d", round($$boxRef{'cooX'} / $gridDeltaX, 1) + 1);
#    my $gridsizeY = sprintf("%d", round($$boxRef{'cooY'} / $gridDeltaY, 1) + 1);
#    my $gridsizeZ = sprintf("%d", round($$boxRef{'cooZ'} / $gridDeltaZ, 1) + 1);
#    my @grid;
#
#    printf("      Create a %dx%dx%d 3D grid: 0%%\r", $gridsizeZ, $gridsizeX, $gridsizeY) if $main::verbose;
#    for (my $z=0; $z<=$gridsizeZ; $z++) {
#        for (my $x=0; $x<=$gridsizeX; $x++) {
#            for (my $y=0; $y<=$gridsizeY; $y++) {
#                $grid[$z][$x][$y]{'type'} = 'VOX';
#            }
#        }
#        printf("      Create a %dx%dx%d 3D grid: %d%%\r", $gridsizeZ, $gridsizeX, $gridsizeY, 100*$z/$gridsizeZ) if $main::verbose;
#    }
#    printf("      Create a %dx%dx%d 3D grid: 100%%\n", $gridsizeZ, $gridsizeX, $gridsizeY) if $main::verbose;
#
#    return (\@grid, $gridsizeZ, $gridsizeX, $gridsizeY);
#}



#sub getMax {
#    return $_[0] if $_[0] > $_[1];
#    return $_[1];
#}
#
#
#
#sub getMin {
#    return $_[0] if $_[0] < $_[1];
#    return $_[1];
#}



#sub analyze {
#    my $atomIdsRef      = shift;
#    my $coordDataRef    = shift;
#    my $suGroupIdsRef   = shift;
#    my $extCavDetection = shift;
#    my $xyProfile       = shift;
#
#    ### Reset grid #############################################################
#    undef(@grid);
##    $gridXMin[$z] = 1000;
##    $gridYMin[$z] = 1000;
##    $gridZMin = 1000;
##    $gridXMax[$z] = 0;
##    $gridYMax[$z] = 0;
##    $gridZMax = 0;
#    ############################################################################
#
#    my $gridVolumeProtein;
#    my $gridVolumeCavity;
#    my $gridVolumeSurface;
#
#
#    ### Initialize grid ########################################################
##    ($gridRef, $gridsizeZ, $gridsizeX, $gridsizeY) = ini3DGrid($boxRef, $gridDelta, 0.1);
#    ############################################################################
#
#
#    ### Create protein grid ####################################################
##    $gridVolumeProtein = protein2Grid($gridRef, $gridsizeZ, $gridsizeX, $gridsizeY, $gridDelta, $atomIdsRef, $coordDataRef);
#    $gridVolumeProtein = protein2Grid($atomIdsRef, $coordDataRef);
##    print "$gridXMin $gridYMin $gridZMin :: $gridXMax $gridYMax $gridZMax\n";
#    ############################################################################
#
#
#    ### Connect subdomains to form a cavity ####################################
##    connectSubunits($gridRef, $gridsizeX, $gridsizeY, $gridDelta, $ndxDataRef, $suGroupIdsRef, $coordsRef, $zMin, $zMax) if $suGroupIdsRef;
#    ############################################################################
#
#
#    ### Detect the cavity ######################################################
#    $gridVolumeCavity = detectCavity($extCavDetection);
#    ############################################################################
#
#
#    ### Detect the surface #####################################################
#    detectSurface();
#    ############################################################################
#
#
#    ### Write out the XY profile ###############################################
##    writeXyProfile ($xyProfile, $gridRef, $gridsizeX, $gridsizeY, $gridDelta) if $xyProfile;
#    ############################################################################
#
#
#    return ($gridVolumeProtein, $gridVolumeCavity, \@grid);
#}



sub detectSurface {
#    print "      Detecting protein surface: 0%\r" if $main::verbose;
#    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
#        for (my $x=$gridXMin; $x<=$gridXMax; $x++) {
#            for (my $y=$gridYMin; $y<=$gridYMax; $y++) {
#                $grid[$z][$x][$y]{'type'} = 'VOX' unless $grid[$z][$x][$y];
#                next unless $grid[$z][$x][$y]{'type'} eq 'PRO';
#                if (getNeighborVox($grid[$z][$x][$y], $grid[$z][$x-1][$y])) {
#                    $gridXMin = $x if $x < $gridXMin;
#                }
#                if (getNeighborVox($grid[$z][$x][$y], $grid[$z][$x+1][$y])) {
#                    $gridXMax = $x if $x > $gridXMax;
#                }
#                if (getNeighborVox($grid[$z][$x][$y], $grid[$z][$x][$y-1])) {
#                    $gridYMin = $y if $y < $gridYMin;
#                }
#                if (getNeighborVox($grid[$z][$x][$y], $grid[$z][$x][$y+1])) {
#                    $gridYMax = $y if $y > $gridYMax;
#                }
#                if (getNeighborVox($grid[$z][$x][$y], $grid[$z][$x+1][$y+1])) {
#                    $gridXMax = $x if $x > $gridXMax;
#                    $gridYMax = $y if $y > $gridYMax;
#                }
#                if (getNeighborVox($grid[$z][$x][$y], $grid[$z][$x+1][$y-1])) {
#                    $gridXMax = $x if $x > $gridXMax;
#                    $gridYMin = $y if $y < $gridYMin;
#                }
#                if (getNeighborVox($grid[$z][$x][$y], $grid[$z][$x-1][$y+1])) {
#                    $gridXMin = $x if $x < $gridXMin;
#                    $gridYMax = $y if $y > $gridYMax;
#                }
#                if (getNeighborVox($grid[$z][$x][$y], $grid[$z][$x-1][$y-1])) {
#                    $gridXMin = $x if $x < $gridXMin;
#                    $gridYMin = $y if $y < $gridYMin;
#                }
#            }
#        }
#        printf("      Detecting protein surface: %d%%      \r", $z/$gridZMax*100) if $main::verbose;
#    }
#    print "      Detecting protein surface: 100%      \n" if $main::verbose;
##    return ($nCavityAreas*$gridDelta2);
}



sub getNeighborVox {
    my $protVoxelRef  = shift;
    my $neighVoxelRef = shift;

    return unless $$neighVoxelRef{'type'};
    return if $$neighVoxelRef{'type'} eq 'PRO';
    return if $$neighVoxelRef{'type'} eq 'INT';

    $$neighVoxelRef{'type'} = 'SUR';
    $$neighVoxelRef{'sum'} += $$protVoxelRef{'sum'};

    return 1;
}


1;
