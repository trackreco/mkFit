#!/usr/bin/perl

@Offs[3] = [ 0, 1, 3, 1, 2, 4, 3, 4, 5 ];
@Offs[4] = [ 0, 1, 3, 6, 1, 2, 4, 7, 3, 4, 5, 8, 6, 7, 8, 9 ];
@Offs[5] = [ 0, 1, 3, 6, 10, 1, 2, 4, 7, 11, 3, 4, 5, 8, 12, 6, 7, 8, 9, 13, 10, 11, 12, 13, 14 ];
@Offs[6] = [ 0, 1, 3, 6, 10, 15, 1, 2, 4, 7, 11, 16, 3, 4, 5, 8, 12, 17, 6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20 ];

$PREF = "      ";
$BR   = " ";
$JOIN = "$BR";
$POST = " ";

$D   = 3;
@Off = @{$Offs[$D]};

$SYMMETRIC = 1;

if ($SYMMETRIC)
{ 

  for (my $i = 0; $i < $D; ++$i)
  {
    for (my $j = 0; $j < $D; ++$j)
    {
      # my $x = $Off[$i * $D + $j];
      my $x = $i * $D + $j;
      print "${PREF}C.fArray[$x * N + n] =${POST}";

      my @sum;

      for (my $k = 0; $k < $D; ++$k)
      {
        my $iko = $Off[$i * $D + $k];
        my $kjo = $Off[$k * $D + $j];

        push @sum, "A.fArray[$iko * N + n] * B.fArray[$kjo * N + n]";
      }
      print join(" +$JOIN", @sum), ";${POST}";
      print "\n";
    }
  }

}
else
{

  for (my $i = 0; $i < $D; ++$i)
  {
    for (my $j = 0; $j < $D; ++$j)
    {
      my $x = $i * $D + $j;
      print "${PREF}C.fArray[$x * N + n] =${POST}";

      my @sum;

      for (my $k = 0; $k < $D; ++$k)
      {
        my $iko = $Off[$i * $D + $k];
        my $kjo = $Off[$k * $D + $j];

        push @sum, "A.fArray[$iko * N + n] * B.fArray[$kjo * N + n]";
      }
      print join(" +$JOIN", @sum), ";${POST}";
      print "\n";
    }
  }

}
