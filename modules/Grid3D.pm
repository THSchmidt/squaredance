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


my @grid;

#my $gridXMin = 1000;
#my $gridXMax = 0;
#
#my $gridYMin = 1000;
#my $gridYMax = 0;

my @gridXMin;
my @gridXMax;
my @gridYMin;
my @gridYMax;
my $gridZMin = 1000;
my $gridZMax = 0;

my $gridDeltaX = 0.5;
my $gridDeltaY = 0.5;
my $gridDeltaZ = 0.5;

my $coordShift = 1; # Shift all coordinates to avoid negative grid points (PBC).


sub getMin {
    my $min = $_[$gridZMin];
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        next unless $_[$z];
        $min = $_[$z] if $_[$z] < $min;
    }
    return $min;
}



sub getMax {
    my $max = $_[$gridZMax];
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        next unless $_[$z];
        $max = $_[$z] if $_[$z] > $max;
    }
    return $max;
}



sub getGridRef {
    return \@grid;
}



sub grid2GroFile {
    my $gridOutFile  = shift;

    my $voxId        = 0;
    my %gridGroData;

    my $gridXMinAbs = getMin(@gridXMin);
    my $gridYMinAbs = getMin(@gridYMin);
    my $gridXMaxAbs = getMax(@gridXMax);
    my $gridYMaxAbs = getMax(@gridYMax);

    for (my $z=$gridZMin; $z<=$gridZMax; $z+=2) {
        for (my $x=$gridXMinAbs; $x<=$gridXMaxAbs; $x+=2) {
            for (my $y=$gridYMinAbs; $y<=$gridYMaxAbs; $y+=2) {
                next unless $grid[$z][$x][$y];
                my $resName = defined($grid[$z][$x][$y]{'sum'}) ? sprintf("%d", $grid[$z][$x][$y]{'sum'}) : $grid[$z][$x][$y]{'type'};
                my %tmpCoords = ('cooX' => $x * $gridDeltaX - $coordShift,
                                 'cooY' => $y * $gridDeltaY - $coordShift,
                                 'cooZ' => $z * $gridDeltaZ - $coordShift);
                push(@{$gridGroData{'atoms'}}, setVoxel(++$voxId, $resName, $grid[$z][$x][$y]{'type'}, $voxId, \%tmpCoords));
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



sub setGridDelta {
    return 0 unless defined $_[0];
    $gridDeltaX = $_[0];
    $gridDeltaY = defined $_[1] ? $_[1] : $_[0];
    $gridDeltaZ = defined $_[2] ? $_[2] : $_[1];
    return 1;
}



sub analyze {
    my $atomIdsRef      = shift;
    my $coordDataRef    = shift;
    my $suGroupIdsRef   = shift;
    my $extCavDetection = shift;
    my $xyProfile       = shift;

    ### Reset grid #############################################################
    undef(@grid);
#    $gridXMin[$z] = 1000;
#    $gridYMin[$z] = 1000;
#    $gridZMin = 1000;
#    $gridXMax[$z] = 0;
#    $gridYMax[$z] = 0;
#    $gridZMax = 0;
    ############################################################################

    my $gridVolumeProtein;
    my $gridVolumeCavity;
    my $gridVolumeSurface;


    ### Initialize grid ########################################################
#    ($gridRef, $gridsizeZ, $gridsizeX, $gridsizeY) = ini3DGrid($boxRef, $gridDelta, 0.1);
    ############################################################################


    ### Create protein grid ####################################################
#    $gridVolumeProtein = protein2Grid($gridRef, $gridsizeZ, $gridsizeX, $gridsizeY, $gridDelta, $atomIdsRef, $coordDataRef);
    $gridVolumeProtein = protein2Grid($atomIdsRef, $coordDataRef);
#    print "$gridXMin $gridYMin $gridZMin :: $gridXMax $gridYMax $gridZMax\n";
    ############################################################################


    ### Connect subdomains to form a cavity ####################################
#    connectSubunits($gridRef, $gridsizeX, $gridsizeY, $gridDelta, $ndxDataRef, $suGroupIdsRef, $coordsRef, $zMin, $zMax) if $suGroupIdsRef;
    ############################################################################


    ### Detect the cavity ######################################################
    $gridVolumeCavity = detectCavity($extCavDetection);
    ############################################################################


    ### Detect the surface #####################################################
    detectSurface();
    ############################################################################


    ### Write out the XY profile ###############################################
#    writeXyProfile ($xyProfile, $gridRef, $gridsizeX, $gridsizeY, $gridDelta) if $xyProfile;
    ############################################################################


    return ($gridVolumeProtein, $gridVolumeCavity, \@grid);
}



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



sub atoms2Grid {
    my $atomsIdsRef  = shift;
    my $coordDataRef = shift;

    my $atomId   = 0;

    print "  ---------------------------------\n  Mapping atoms to the grid\r";

    foreach (@{$atomsIdsRef}) {
        next unless $$coordDataRef[$_]{'cooZ'};

        my $element = substr($$coordDataRef[$_]{'atomName'}, 0, 1);
#        my $element = $$coordDataRef[$_]{'atomName'};
        my $radius  = PTE::getRadius($element);
        my $radius2 = $radius * $radius;


        my $tmpGridX = sprintf("%d", round(($$coordDataRef[$_]{'cooX'} + $coordShift) / $gridDeltaX, 1));
        my $tmpGridY = sprintf("%d", round(($$coordDataRef[$_]{'cooY'} + $coordShift) / $gridDeltaY, 1));
        my $tmpGridZ = sprintf("%d", round(($$coordDataRef[$_]{'cooZ'} + $coordShift) / $gridDeltaZ, 1));
        my $subrangeX = $radius > $gridDeltaX ? sprintf("%d", round($radius / $gridDeltaX, 1)) : 0;
        my $subrangeY = $radius > $gridDeltaY ? sprintf("%d", round($radius / $gridDeltaY, 1)) : 0;
        my $subrangeZ = $radius > $gridDeltaZ ? sprintf("%d", round($radius / $gridDeltaZ, 1)) : 0;

        for (my $z=($tmpGridZ-$subrangeZ); $z<=($tmpGridZ+$subrangeZ); $z++) {
            for (my $x=($tmpGridX-$subrangeX); $x<=($tmpGridX+$subrangeX); $x++) {
                for (my $y=($tmpGridY-$subrangeY); $y<=($tmpGridY+$subrangeY); $y++) {
                    my $dx = $$coordDataRef[$_]{'cooX'} + $coordShift - $x * $gridDeltaX;
                    my $dy = $$coordDataRef[$_]{'cooY'} + $coordShift - $y * $gridDeltaY;
                    my $dz = $$coordDataRef[$_]{'cooZ'} + $coordShift - $z * $gridDeltaZ;
                    next if ($dx*$dx + $dy*$dy) > $radius2;
                    next if ($dx*$dx + $dz*$dz) > $radius2;
                    next if ($dy*$dy + $dz*$dz) > $radius2;

                    $grid[$z][$x][$y]{'type'} = 'ATM';
#                    $grid[$z][$x][$y]{'type'} = $$coordDataRef[$_]{'resName'};

                    $gridXMin[$z] = $x if !defined $gridXMin[$z] || $x < $gridXMin[$z];
                    $gridXMax[$z] = $x if !defined $gridXMax[$z] || $x > $gridXMax[$z];
                    $gridYMin[$z] = $y if !defined $gridYMin[$z] || $y < $gridYMin[$z];
                    $gridYMax[$z] = $y if !defined $gridYMax[$z] || $y > $gridYMax[$z];
                    $gridZMin = $z if $z < $gridZMin;
                    $gridZMax = $z if $z > $gridZMax;
                }
            }
        }
        printf("  Mapping atoms to the grid: %d%%\r", ++$atomId*100/@{$atomsIdsRef}) if $main::verbose;
    }

    printf("  Mapping atoms to the grid: Finished\n  ---------------------------------\n\n") if $main::verbose;

    return 1;
}



sub extendGridXY {
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        $gridXMin[$z]--;
        $gridYMin[$z]--;
        $gridXMax[$z]++;
        $gridYMax[$z]++;
    }
}
    ### Plus imaginary space for surface voxels and cavity analysis ############
    ############################################################################



sub getSas {
    my $probeRadius = shift;

    return if $probeRadius == 0;
    $probeRadius = 0.14 unless defined $probeRadius;

    
}


sub detectCavity {
    my $maxRadius = shift;

    $maxRadius = 0.2 unless defined $maxRadius; # This could also be 0 to washing out cavities with even small gaps to the exterior.
    my $maxRadiusGridX = sprintf("%d", round($maxRadius / $gridDeltaX, 1));
    my $maxRadiusGridY = sprintf("%d", round($maxRadius / $gridDeltaY, 1));

    my $nCavityAreas = 0;
    my $gridDelta2   = $gridDeltaX * $gridDeltaY;

    extendGridXY();

    print "  ---------------------------------\n  Detecting internal cavities\r" if $main::verbose;

    ### Invert existing grid (define void voxels as cavity) ####################
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        for (my $x=$gridXMin[$z]; $x<=$gridXMax[$z]; $x++) {
            for (my $y=$gridYMin[$z]; $y<=$gridYMax[$z]; $y++) {
                next if defined $grid[$z][$x][$y] && $grid[$z][$x][$y];
                $grid[$z][$x][$y]{'type'} = 'INT';
                $nCavityAreas++;
            }
            printf("  Detecting internal cavities: %.4f nm^2 possible      \r", $nCavityAreas*$gridDelta2) if $main::verbose;
        }
    }
    ############################################################################


    ### Washing out the cavities ###############################################
    for (my $z=$gridZMin; $z<=$gridZMax; $z++) {
        for (my $i=0; $i<2; $i++) {
            my $foundExcl1 = 1;
            my $foundExcl2 = 1;
            for (my $x=$gridXMin[$z]; $x<=$gridXMax[$z]; $x++) { # W -> E.
                $nCavityAreas -= $foundExcl1 = washingOutY($x, $z, $x-1, $grid[$z], $maxRadiusGridX, $maxRadiusGridY);
                printf("  Detecting internal cavities: %.4f nm^2 possible      \r", $nCavityAreas*$gridDelta2) if $main::verbose;
            }

            for (my $x=$gridXMax[$z]; $x>=$gridXMin[$z]; $x--) { # E -> W.
                $nCavityAreas -= $foundExcl2 = washingOutY($x, $z, $x+1, $grid[$z], $maxRadiusGridX, $maxRadiusGridY);
                printf("  Detecting internal cavities: %.4f nm^2 possible      \r", $nCavityAreas*$gridDelta2) if $main::verbose;
            }
            $i = 0 if $foundExcl1 || $foundExcl2;
        }
    }
    ############################################################################

    printf("  Detecting internal cavities: Finished\n  ---------------------------------\n\n", $nCavityAreas*$gridDelta2) if $main::verbose;

    return ($nCavityAreas*$gridDelta2);
}



sub round {
    my ( $num, $prec ) = @_;
    return int( $num / $prec + 0.5 - ( $num < 0 ) ) * $prec;
}



sub searchOverlaps {
    my $gridSliceRef   = shift;
    my $x              = shift;
    my $y              = shift;
    my $maxRadiusGridX = shift;
    my $maxRadiusGridY = shift;

    return 'VOX' unless $maxRadiusGridX && $maxRadiusGridY;
    return 'ASA' if $$gridSliceRef[$x+$maxRadiusGridX][$y] && $$gridSliceRef[$x+$maxRadiusGridX][$y]{'type'} eq 'ATM';
    return 'ASA' if $$gridSliceRef[$x][$y+$maxRadiusGridY] && $$gridSliceRef[$x][$y+$maxRadiusGridY]{'type'} eq 'ATM';
    return 'ASA' if $$gridSliceRef[$x-$maxRadiusGridX][$y] && $$gridSliceRef[$x-$maxRadiusGridX][$y]{'type'} eq 'ATM';
    return 'ASA' if $$gridSliceRef[$x][$y-$maxRadiusGridY] && $$gridSliceRef[$x][$y-$maxRadiusGridY]{'type'} eq 'ATM';
    return 'VOX';
}



sub washingOutY {
    my $x            = shift;
    my $z            = shift;
    my $neighbCellX  = shift;
    my $gridSliceRef = shift;
    my $maxRadiusGridX = shift;
    my $maxRadiusGridY = shift;

    my $delAreas        = 0;

    for (my $y=$gridYMin[$z]; $y<=$gridYMax[$z]; $y++) { # S -> N.
        my $neighbCellY = $y - 1;
        next unless $$gridSliceRef[$x][$y]{'type'} eq 'INT'; # Next if this grid point is not a possible cavity.
        if (!$$gridSliceRef[$x][$neighbCellY] || $$gridSliceRef[$x][$neighbCellY]{'type'} eq 'VOX') { # If the neighbored cell is not protein or not a possible cavity...
#            $$gridSliceRef[$x][$y]{'type'} = 'VOX'; # ...exclude this from a possible cavity.
            $$gridSliceRef[$x][$y]{'type'} = searchOverlaps($gridSliceRef, $x, $y, $maxRadiusGridX, $maxRadiusGridY); # ...exclude this from a possible cavity.
            if ($$gridSliceRef[$x][$y]{'type'} eq 'VOX') {
                $delAreas++;
                next;
            }
        }
        if (!$$gridSliceRef[$neighbCellX][$y] || $$gridSliceRef[$neighbCellX][$y]{'type'} eq 'VOX') {
#            $$gridSliceRef[$x][$y]{'type'} = 'VOX';
            $$gridSliceRef[$x][$y]{'type'} = searchOverlaps($gridSliceRef, $x, $y, $maxRadiusGridX, $maxRadiusGridY); # ...exclude this from a possible cavity.
            if ($$gridSliceRef[$x][$y]{'type'} eq 'VOX') {
                $delAreas++;
                next;
            }
        }
    }

    for (my $y=$gridYMax[$z]; $y>=$gridYMin[$z]; $y--) { # N -> S.
        my $neighbCellY = $y + 1;
        next unless $$gridSliceRef[$x][$y]{'type'} eq 'INT';
        if (!$$gridSliceRef[$x][$neighbCellY] || $$gridSliceRef[$x][$neighbCellY]{'type'} eq 'VOX') {
#            $$gridSliceRef[$x][$y]{'type'} = 'VOX';
            $$gridSliceRef[$x][$y]{'type'} = searchOverlaps($gridSliceRef, $x, $y, $maxRadiusGridX, $maxRadiusGridY); # ...exclude this from a possible cavity.
            if ($$gridSliceRef[$x][$y]{'type'} eq 'VOX') {
                $delAreas++;
                next;
            }
        }
        if (!$$gridSliceRef[$neighbCellX][$y] || $$gridSliceRef[$neighbCellX][$y]{'type'} eq 'VOX') {
#            $$gridSliceRef[$x][$y]{'type'} = 'VOX';
            $$gridSliceRef[$x][$y]{'type'} = searchOverlaps($gridSliceRef, $x, $y, $maxRadiusGridX, $maxRadiusGridY); # ...exclude this from a possible cavity.
            if ($$gridSliceRef[$x][$y]{'type'} eq 'VOX') {
                $delAreas++;
                next;
            }
        }
    }

    return $delAreas;
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



sub connectSubunits {
    my $gridRef       = shift;
    my $gridsizeX     = shift;
    my $gridsizeY     = shift;
    my $gridDelta     = shift;
    my $ndxDataRef    = shift;
    my $suGroupIdsRef = shift;
    my $coordsRef     = shift;
    my $zMin          = shift;
    my $zMax          = shift;

    my @gridSuGeoCenter;


    print "      Detecting subunit connections: 0.0000 nm^2 possible\r" if $main::verbose;

    ### Calculate center of mass of each protein subunit in that slice #########
    foreach my $suGroupId (@{$suGroupIdsRef}) {
        my %tmpSum = ('cooX' => 0, 'cooY' => 0);
        my $nAtoms = 0;
        foreach (@{$$ndxDataRef[$suGroupId]{'atoms'}}) {
            next unless $$coordsRef[$_]{'cooZ'};
            next if $$coordsRef[$_]{'cooZ'} > $zMax;
            next if $$coordsRef[$_]{'cooZ'} < $zMin;
            $tmpSum{'cooX'} += $$coordsRef[$_]{'cooX'};
            $tmpSum{'cooY'} += $$coordsRef[$_]{'cooY'};
            $nAtoms++;
        }
        next unless $nAtoms;
        $tmpSum{'cooX'} = int(($tmpSum{'cooX'}/$nAtoms)/$gridDelta);
        $tmpSum{'cooY'} = int(($tmpSum{'cooY'}/$nAtoms)/$gridDelta);
        push(@gridSuGeoCenter, \%tmpSum);
    }
    ############################################################################


    ### Print the geometrical center of each subunit ###########################
    for (my $i=0; $i<@gridSuGeoCenter; $i++) {
        getGridSlope($gridRef, $gridSuGeoCenter[$i-1], $gridSuGeoCenter[$i]);
    }
    ############################################################################
}



sub getGridSlope {
    my $gridRef   = shift;
    my $vecARef   = shift;
    my $vecBRef   = shift;

    my $slope = ($$vecBRef{'cooY'} - $$vecARef{'cooY'}) / ($$vecBRef{'cooX'} - $$vecARef{'cooX'});
    my $xMin = $$vecARef{'cooX'} < $$vecBRef{'cooX'} ? $$vecARef{'cooX'} : $$vecBRef{'cooX'};
    my $xMax = $$vecARef{'cooX'} > $$vecBRef{'cooX'} ? $$vecARef{'cooX'} : $$vecBRef{'cooX'};

    my $y = $$vecARef{'cooY'};
    for (my $x=$xMin; $x<=$xMax; $x++) {
        $y += $slope;
        my $intY = int($y);
        next if $$gridRef[$x][$intY];
        $$gridRef[$x][$intY] = 3;

        $$gridRef[$x+1][$intY] = 3;
        $$gridRef[$x-1][$intY] = 3;
        $$gridRef[$x][$intY+1] = 3;
        $$gridRef[$x][$intY-1] = 3;
        $$gridRef[$x+1][$intY+1] = 3;
        $$gridRef[$x+1][$intY-1] = 3;
        $$gridRef[$x-1][$intY+1] = 3;
        $$gridRef[$x-1][$intY-1] = 3;
    }
}



sub ini3DGrid {
    my $boxRef     = shift;
    my $gridDeltaZ = shift;
    my $gridDeltaX = shift;
    my $gridDeltaY = shift;

    $gridDeltaZ = 0.5 unless $gridDeltaZ;
    $gridDeltaX = $gridDeltaZ unless $gridDeltaX;
    $gridDeltaY = $gridDeltaX unless $gridDeltaY;

    my $gridsizeX = sprintf("%d", round($$boxRef{'cooX'} / $gridDeltaX, 1) + 1);
    my $gridsizeY = sprintf("%d", round($$boxRef{'cooY'} / $gridDeltaY, 1) + 1);
    my $gridsizeZ = sprintf("%d", round($$boxRef{'cooZ'} / $gridDeltaZ, 1) + 1);
    my @grid;

    printf("      Create a %dx%dx%d 3D grid: 0%%\r", $gridsizeZ, $gridsizeX, $gridsizeY) if $main::verbose;
    for (my $z=0; $z<=$gridsizeZ; $z++) {
        for (my $x=0; $x<=$gridsizeX; $x++) {
            for (my $y=0; $y<=$gridsizeY; $y++) {
                $grid[$z][$x][$y]{'type'} = 'VOX';
            }
        }
        printf("      Create a %dx%dx%d 3D grid: %d%%\r", $gridsizeZ, $gridsizeX, $gridsizeY, 100*$z/$gridsizeZ) if $main::verbose;
    }
    printf("      Create a %dx%dx%d 3D grid: 100%%\n", $gridsizeZ, $gridsizeX, $gridsizeY) if $main::verbose;

    return (\@grid, $gridsizeZ, $gridsizeX, $gridsizeY);
}



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


1;
