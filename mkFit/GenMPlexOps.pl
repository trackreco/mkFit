#!/usr/bin/perl

use lib "../Matriplex";

use GenMul;
use warnings;

my $DIM = 6;

### Propagate Helix To R -- final similarity, two ops.

# outErr = errProp * outErr * errPropT
#   outErr is symmetric


$errProp = new GenMul::Matrix('name'=>'a', 'M'=>$DIM, 'N'=>$DIM);
$errProp->set_pattern(<<"FNORD");
x x 0 x x 0
x x 0 x x 0
x x 1 x x x
x x 0 x x 0
x x 0 x x 0
0 0 0 0 0 1
FNORD

$outErr = new GenMul::MatrixSym('name'=>'b', 'M'=>$DIM, 'N'=>$DIM);

$temp   = new GenMul::Matrix('name'=>'c', 'M'=>$DIM, 'N'=>$DIM);


$errPropT = new GenMul::MatrixTranspose($errProp);
$errPropT->print_info();
$errPropT->print_pattern();

# ----------------------------------------------------------------------

$m = new GenMul::Multiply;

# outErr and c are just templates ...

$m->dump_multiply_std_and_intrinsic("MultHelixProp.ah",
                                    $errProp, $outErr, $temp);

$temp  ->{name} = 'b';
$outErr->{name} = 'c';

### XXX fix this ... in accordance with what is in Propagation.cc
$m->dump_multiply_std_and_intrinsic("MultHelixPropTransp.ah",
                                    $temp, $errPropT, $outErr);



##############################
### updateParameters       ###
##############################

#declared first on its own because propErr sees many uses
my $propErr_M = 6;
$propErr = new GenMul::MatrixSym('name' => 'a',
                                 'M'    => $propErr_M); #will have to remember to re'name' it based on location in function

my $propErrT_M = 6;
$propErrT = new GenMul::MatrixTranspose($propErr); #will have to remember to re'name' it based on location in function



### kalmanGain =  = propErr * (projMatrixT * resErrInv)
$resErrInv = new GenMul::MatrixSym('name'=>'b', 'M'=>3, 'N'=>3);

$kalmanGain = new GenMul::Matrix('name'=>'c', 'M' => 6, 'N' => 3);

{
  my $m_kg = new GenMul::Multiply('no_size_check' => 1);

  $m_kg->dump_multiply_std_and_intrinsic("upParam_MultKalmanGain.ah",
                                         $propErr, $resErrInv, $kalmanGain);
}


### updatedErrs = propErr - propErr^T * simil * propErr
# Going to skip the subtraction for now
my $simil_M = 6;
$simil = new GenMul::MatrixSym('name'=>'a', 'M'=>$simil_M);
$simil->set_pattern(<<"FNORD");
x
x x
x x x
0 0 0 0
0 0 0 0 0
0 0 0 0 0 0
FNORD

$propErr->{name} = 'b';

my $temp_simil_x_propErr_M = 6;
my $temp_simil_x_propErr_N = 6;
$temp_simil_x_propErr = new GenMul::Matrix('name'=>'c',
                                           'M'=>$temp_simil_x_propErr_M,
                                           'N'=>$temp_simil_x_propErr_N);

$m->dump_multiply_std_and_intrinsic("upParam_simil_x_propErr.ah",
                                    $simil, $propErr, $temp_simil_x_propErr);

$temp_simil_x_propErr->{name} = 'b';									 
$temp_simil_x_propErr->set_pattern(<<"FNORD");
x x x x x x
x x x x x x
x x x x x x
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
FNORD

#? This one is symmetric but the output can't handle it... need to fix
#$temp_propErrT_x_simil_propErr = new GenMul::MatrixSym('name'=>'c', 'M'=>$propErrT_M, 'N'=>$temp_simil_x_propErr_N);


$temp_propErrT_x_simil_propErr = new GenMul::MatrixSym('name'=>'c', 'M'=>$propErrT_M);

$m->dump_multiply_std_and_intrinsic("upParam_propErrT_x_simil_propErr.ah",
                                    $propErrT, $temp_simil_x_propErr, $temp_propErrT_x_simil_propErr);
									

{
  my $temp = new GenMul::MatrixSym('name' => 'c', 'M' => 6);

  my $m_kg = new GenMul::Multiply('no_size_check' => 1);

  $kalmanGain->{name} = 'a';

  $m_kg->dump_multiply_std_and_intrinsic("upParam_kalmanGain_x_propErr.ah",
                                         $kalmanGain, $propErr, $temp);
}
