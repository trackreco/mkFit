#!/usr/bin/perl

=head1

./chichi.pl std
./chichi.pl wfix --kludge-cms-hit-errors

=cut

@files = map { "/bar/mic/mu_$_-1000-10.bin-5" } qw(brl ecn ecp trn trp);

$flds_hit = "layer/I:chi2/F:x_h:y_h:z_h:r_h:ex_h:ey_h:ez_h:x_t:y_t:z_t:r_t:ex_t:ey_t:ez_t:pt:phi:theta:phi_h:phi_t:ephi_h:ephi_t";
$flds_trk = "n_hits/I:chi2/F:chi2pdof:pt:phi:theta";

$base_opts   = "--build-ce --geom CMS-2017 --num-events 990 --backward-fit";

#-----------------------------------------------------------------------

die "Usage: $0 output-stem [additional mkFit arguments]" if ($#ARGV < 0);

$stem       = shift @ARGV;
$extra_opts = join(" ", @ARGV);

print "$0 using stem='$stem', extra_opts='$extra_opts'\n";

#-----------------------------------------------------------------------

sub make_with_opts
{
  my $name = shift;
  my $opts = shift;

  my $hit = "$stem-hit-$name.rtt", $trk = "$stem-trk-$name.rtt";

  unlink $hit, $trk;

  open H, ">$hit"; print H $flds_hit, "\n";
  open T, ">$trk"; print T $flds_trk, "\n";

  print "make_with_opts $name - opts='$opts'\n";

  for $f (@files)
  {
    print "    processing $f\n";

    open MKF, "./mkFit --read --file-name $f ${base_opts} ${extra_opts} ${opts} |";

    while ($_ = <MKF>)
    {
      print H if s/^CHIHIT //o;
      print T if s/^CHITRK //o;
    }

    close MKF;
  }

  close H;
  close T;
}

#-----------------------------------------------------------------------

make_with_opts("pure", "--cmssw-seeds");

make_with_opts("sim", "");
