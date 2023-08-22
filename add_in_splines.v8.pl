#!/usr/bin/perl -w

use strict;
use Getopt::Long;

## Goal: Recompute spline points for all SplineRoi elements. 
## Discard existing spline points.
## Created: March 23, 2017 by TS
## NOTE: Have to add in the first polygon point at the end before calling poly_to_spline.
## Version 2 - April 17, 2017 by TS
## ** adding in support for spline points in PolySplines (open isopters)
## Version 2 - April 21, 2017 by TS
## ** add in support for specifying number of spline points
## **  as total number --spline-points, or relative to # of polygon points as --increase-points

##                       <name>I3e-open                        </name>
##                       <value>I3e-open                        </value>
##                       <controltype>PolysplineRoi</controltype>
##                       <roi>
##                                <polyline>
##                                        <npoints>35</npoints>
##                                        <xpoints>
##                                                <int>449</int>



my $NUM_SPLINE_POINTS = 750;
my $file = "";
my $SCALE_FROM_POLYGON = "";

GetOptions("file=s" => \$file, "spline-points=i" => \$NUM_SPLINE_POINTS, "increase-points=f" => \$SCALE_FROM_POLYGON) or die("Failed to parse arguments");

if($file eq "") {
  print STDERR "ERROR: You must specify a file\n";
  Usage();
  exit();
}

#print STDERR "DEBUG: SCALE_FROM_POLYGON = $SCALE_FROM_POLYGON\n";


#my $in_annotation = 0;
my $in_spline_roi = 0;
my $in_polygon_roi = 0;
my $in_polygon_x = 0;
my $in_polygon_y = 0;
my $in_spline_x = 0;
my $in_spline_y = 0;
#my $roi_type = "";
my @x_poly = ();
my @y_poly = ();
my @x_spline = ();
my @y_spline = ();
my $done_for_this_isopter = 0;
my $length = 0;
my @length = ();
my $current_point_x = 0;
my $previous_point_x = 0;
my $current_point_y = 0;
my $previous_point_y = 0;

## open file
open(IN,$file) or die("Failed to open $file");
while(my $line = <IN>) {

  if($line =~ m/<spline>/) {
    $in_spline_roi = 1;
    next();
  }
  if($line =~ m/<\/polygon>/) {
    print $line;
    ## process polygon to spline
    push @x_poly, $x_poly[0];
    push @y_poly, $y_poly[0];
    $length = @x_poly;
    @length = (1..($length-1));
    print STDERR "X_POLY: @x_poly\n";
    #print "@length\n";
    for (my $i=1; $i< $length; $i++) {
        #print "@x_poly[($i)]\n";
        $current_point_x = $x_poly[$i];
        $previous_point_x = $x_poly[($i-1)];
        $current_point_y = $y_poly[$i];
        $previous_point_y = $y_poly[($i-1)];
        #print "current_point_x: $current_point_x\n";
        if ($current_point_x == $previous_point_x){
            if ($current_point_y == $previous_point_y){
                splice  @x_poly, $i, 1;
                splice  @y_poly, $i, 1;
		$length--;
            }
        }
        #print "$i\n";
    }
    #print "length: $length\n";
    #print "x_poly\n";
    print STDERR "X_POLY: @x_poly\n";
    poly_to_spline(\@x_poly, \@y_poly, 0);

    ## re-set variables
    $in_polygon_roi = $in_polygon_x = $in_polygon_y = 0;
    @x_poly = ();
    @y_poly = ();
    $in_spline_roi = 0;
    $done_for_this_isopter = 1;
  } elsif($line =~ m/<\/polyline>/) {
    print STDERR "Inside matched </polyline>\n";
    print $line;
    ## process polyline to spline
    poly_to_spline(\@x_poly, \@y_poly, 1);

    ## re-set variables
    $in_polygon_roi = $in_polygon_x = $in_polygon_y = 0;
    @x_poly = ();
    @y_poly = ();
    $in_spline_roi = 0;
  } else {
    if($in_spline_roi == 0) {
      print $line;
    }

    if($in_spline_roi == 1) {

    }

    if($in_polygon_x == 1) {
      if($line =~ m/<int>(.*?)<\/int>/) {
        push @x_poly, $1;
      }
    } elsif($in_polygon_y == 1) {
      if($line =~ m/<int>(.*?)<\/int>/) {
        push @y_poly, $1;
      }
    }

    if($in_spline_x == 1) {
      if($line =~ m/<int>(.*?)<\/int>/) {
        push @x_spline, $1;
      }
    } elsif($in_spline_y == 1) {
      if($line =~ m/<int>(.*?)<\/int>/) {
        push @y_spline, $1;
      }
    }
  }
  if($line =~ m/<\/spline>/) {
    $in_spline_roi = 0;

    ## if @X_poly == emtpy (no polygon points) print out spline points as both polygon and spline points
    if($#x_poly < 0 && $done_for_this_isopter == 0) {
      reprint_spline_and_copy_as_poly(\@x_spline, \@y_spline);
    }

    @x_spline = ();
    @y_spline = ();

    $done_for_this_isopter = 0;
  }

  if($line =~ m/<polygon>|<polyline>/) {
    $in_polygon_roi = 1;
  } elsif($line =~ m/<xpoints>/) {
    if($in_polygon_roi == 1) {
      $in_polygon_x = 1;
    }
    if($in_spline_roi == 1) {
      $in_spline_x = 1;
    }
  } elsif($line =~ m/<\/xpoints>/) {
    if($in_polygon_roi == 1) {
      $in_polygon_x = 0;
    }
    if($in_spline_roi == 1) {
      $in_spline_x = 0;
    }
  } elsif($line =~ m/<ypoints>/) {
    if($in_polygon_roi == 1) {
      $in_polygon_y = 1;
    }
    if($in_spline_roi == 1) {
      $in_spline_y = 1;
    }
  } elsif($line =~ m/<\/ypoints>/) {
    if($in_polygon_roi == 1) {
      $in_polygon_y = 0;
    }
    if($in_spline_roi == 1) {
      $in_spline_y = 0;
    }
  } elsif($line =~ m/<\/polygon>/) {
    $in_polygon_roi = 0;
  }



}
close(IN);

exit();

## DEBUG case
#my @x_poly = (414, 409, 395, 389, 389, 402, 414);
#my @y_poly = ( 331, 314, 316, 328, 341, 352, 331);
#print "X = @x_poly\n";
#print "Y = @y_poly\n";
#poly_to_spline(\@x_poly, \@y_poly);


exit();

##############################################################
#### nNodes = # of points/corners in the polygon
#### splinePoints = # of points at which the polygon should be evaluated (typ. 750)

## if it is a polyline, then need to add the final polygon points to the end

sub poly_to_spline {
  my $x_poly_ref = shift;
  my $y_poly_ref = shift;
  my $is_polyline = 0;

  #print STDERR "DEBUG: Inside poly_to_spline\n";
  #print "x_poly_ref";
  #print "$x_poly_ref";
  #print "y_poly_ref";
  #print "$y_poly_ref";


  my $nNodes = @{$x_poly_ref};
  my $splinePoints = $NUM_SPLINE_POINTS;
  if($SCALE_FROM_POLYGON ne "") {
    #$splinePoints = int($SCALE_FROM_POLYGON * ($nNodes-1));
    $splinePoints = int(($SCALE_FROM_POLYGON * ($nNodes-1)) + $nNodes);
  }

  my @index_arr = ();
  for(my $ii = 0; $ii < $nNodes; $ii++) {
    $index_arr[$ii] = $ii;
  }

  #### Allocate set of X and Y spline points
  my @x_spline_points = ();
  my @y_spline_points = ();

  #### Allocate SplineFitter for X and Y 
  my @splineFitterX = ();
  my @splineFitterY = ();
  init_spline_fitter(\@splineFitterX, \@index_arr, $x_poly_ref, $nNodes);
  init_spline_fitter(\@splineFitterY, \@index_arr, $y_poly_ref, $nNodes);

  #### Evaluate the spline at #splinePoints points
  ## ORIGINAL
  my $scale = (1.0 * ($nNodes-1)) / ($splinePoints-1);
  ## testing TODD
  #my $scale = (1.0 * ($nNodes)) / ($splinePoints-1);
  for(my $i=0; $i<$splinePoints; $i++) {
    my $xvalue = 1.0 * $i * $scale;
    my $xpoint = int(eval_spline(\@splineFitterX, \@index_arr, $x_poly_ref, $nNodes, $xvalue));
    my $ypoint = int(eval_spline(\@splineFitterY, \@index_arr, $y_poly_ref, $nNodes, $xvalue));
    ## add to growing spline point array
    push @x_spline_points, $xpoint;
    push @y_spline_points, $ypoint;
  }


  ## add on terminal point (test Todd for open polylines - 20201118)
  if($is_polyline == 1) {
    push @x_spline_points, $x_poly_ref->[$#{$x_poly_ref}];
    push @y_spline_points, $y_poly_ref->[$#{$y_poly_ref}];
  }

  #print STDERR "X: @x_spline_points\n";
  #print STDERR "Y: @y_spline_points\n";

  #### PrintOut the Spline
  print "\t\t\t\t<spline>\n";
  print "\t\t\t\t\t<npoints>$splinePoints</npoints>\n";
  print "\t\t\t\t\t<xpoints>\n";
  for(my $ii=0; $ii <= $#x_spline_points - 1; $ii++) {
    print "\t\t\t\t\t\t<int>$x_spline_points[$ii]</int>\n";
  }
  if($is_polyline == 1) {
    print "\t\t\t\t\t\t<int>$x_spline_points[$#x_spline_points]</int>\n";
  }
  print "\t\t\t\t\t</xpoints>\n";
  print "\t\t\t\t\t<ypoints>\n";
  for(my $ii=0; $ii <= $#y_spline_points - 1; $ii++) {
    print "\t\t\t\t\t\t<int>$y_spline_points[$ii]</int>\n";
  }
  if($is_polyline == 1) {
    print "\t\t\t\t\t\t<int>$y_spline_points[$#y_spline_points]</int>\n";
  }
  print "\t\t\t\t\t</ypoints>\n";
  print "\t\t\t\t</spline>\n";
}

sub reprint_spline_and_copy_as_poly {
  my $xarr_ref = shift;
  my $yarr_ref = shift;

  my $arr_len = @{$xarr_ref};
  ## print Polygon
  print "\t\t\t\t<polygon>\n";
  print "\t\t\t\t\t<npoints>$arr_len</npoints>\n";
  print "\t\t\t\t\t<xpoints>\n";
  for(my $ii=0; $ii<$arr_len; $ii++) {
    print "\t\t\t\t\t\t<int>$xarr_ref->[$ii]</int>\n";
  }
  print "\t\t\t\t\t</xpoints>\n";
  print "\t\t\t\t\t<ypoints>\n";
  for(my $ii=0; $ii<$arr_len; $ii++) {
    print "\t\t\t\t\t\t<int>$yarr_ref->[$ii]</int>\n";
  }
  print "\t\t\t\t\t</ypoints>\n";
  print "\t\t\t\t</polygon>\n";

  ## print Spline
  print "\t\t\t\t<spline>\n";
  print "\t\t\t\t\t<npoints>$arr_len</npoints>\n";
  print "\t\t\t\t\t<xpoints>\n";
  for(my $ii=0; $ii<$arr_len; $ii++) {
    print "\t\t\t\t\t\t<int>$xarr_ref->[$ii]</int>\n";
  }
  print "\t\t\t\t\t</xpoints>\n";
  print "\t\t\t\t\t<ypoints>\n";
  for(my $ii=0; $ii<$arr_len; $ii++) {
    print "\t\t\t\t\t\t<int>$yarr_ref->[$ii]</int>\n";
  }
  print "\t\t\t\t\t</ypoints>\n";
  print "\t\t\t\t</spline>\n";



}

##### Spline Fitter pieces

##init_spline_fitter(\%splineFitterY, \@index_arr, $y_poly_ref, $nNodes);

sub init_spline_fitter {
  my $fitter_ref = shift;   # array ref to store y2
  my $x = shift;
  my $y = shift;
  my $n = shift;

  #print("DEBUG: Inside init_spline_fitter\n");

  #my $i = my $k = 0;
  my $p = my $qn = my $sig = my $un = 0.0;
  my @y2 = ();
  my @u = ();
  for(my $i = 0; $i < $n; $i++) {
    push @y2, 0;
    push @u, 0;
  }

  for(my $i = 1; $i < $n-1; $i++) {
    my $sig = (1.0 * ($x->[$i] - $x->[$i-1])) / (1.0*($x->[$i+1] - $x->[$i-1]));
    my $p = ($sig * $y2[$i-1]) + 2.0;
    $y2[$i] = ($sig-1.0) / $p;
    $u[$i] = ((1.0*($y->[$i+1] - $y->[$i])) / ($x->[$i+1] - $x->[$i])) 
              - ((1.0*($y->[$i] - $y->[$i-1])) / ($x->[$i] - $x->[$i-1]));
    $u[$i] = (6.0 * $u[$i] / ($x->[$i+1] - $x->[$i-1]) - $sig*$u[$i-1])/$p;
  }
  $y2[$n-1] = ($un-$qn*$u[$n-2]) / ($qn*$y2[$n-2]+1.0);

  for (my $k=$n-2; $k>=0; $k--) {
    $y2[$k] = $y2[$k]*$y2[$k+1]+$u[$k];
  }

  for(my $ii=0; $ii<= $#y2; $ii++) {
    push @{$fitter_ref}, $y2[$ii];
  }
  #$fitter_ref->{Y2} = (@y2);
}

## $xpoint = eval_spline(\%splineFitterX, \@index_arr, \@x_points, $nNodes,$xvalue);
sub eval_spline {
  my $spline_fitter_ref = shift;  ## array ref to y2 data
  my $x = shift;
  my $y = shift;
  my $n = shift;
  my $xp = shift;

  #print("DEBUG: Inside eval_spline, xvalue=$xp\n");
  my $klo = 0;
  my $khi = $n-1;
  my $k = 0;
  my $y2_ref = $spline_fitter_ref;
  #print "DEBUG: Y2_ref = $y2_ref\n";

  while($khi-$klo > 1) {
    $k = ($khi+$klo) >> 1;
    if($x->[$k] > $xp) {
      $khi = $k;
    } else {
      $klo=$k
    }
  }

  my $h = $x->[$khi] - $x->[$klo];
  if($h == 0.0) {
    return(0.0);
  }
  my $a = ($x->[$khi] - $xp)/$h;
  my $b = ($xp - $x->[$klo])/$h;

  #if($y2 == NULL) {
  #  return(0.0);
  #}

  my $u = $a * $y->[$klo];
  my $v = $b * $y->[$khi];
  my $y2klo = $y2_ref->[$klo];
  my $y2khi = $y2_ref->[$khi];

  my $w = (($a*$a*$a-$a)*$y2klo + ($b*$b*$b-$b)*$y2khi)*($h*$h) / 6.0;
  my $r = $u + $v + $w;

  return($r);
}

sub Usage {
  print STDERR "$0 --file <FILE> [--spline-points N | --increase-points F]\n";  
}
