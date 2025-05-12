$Global

$Prob
- Busuflan final model - 2 Compartment model with proportional error

$CMT  @annotated
CENT : Central compartment (mg)
PERI : Peripheral (mg)

$PARAM @annotated
CL        :  0.177   : Clearance (volume/time)
VC        :  24.92   : Central volume (volume)
VP        :  14.79   : Peripheral volume (volume)
Q         :  0.344   : Intercompartmental clearance (volume/time)
covBSAVC  :  1.49    : BSA on VC
covBSACL  :  1.03    : BSA on CL
covSEXVP  :  0.757   : Sex on VP
Sex       :  0       : Exponent of weight effect on CL
BSA       :  2.01    : Median BSA

$MAIN
double CLCOV = CL*pow(BSA/2.01, covBSACL);
double CLi   = CLCOV * exp(ECL);
double VCCOV = (VC * exp(EVC))*pow(BSA/2.01, covBSAVC);
double QCOV  = Q;
double VPCOV = (VP * exp(EVP))*pow(covSEXVP, Sex);
double K20   = CLCOV/VCCOV;
double K23   = QCOV/VCCOV;
double K32   = QCOV/VPCOV;
 
$ODE
dxdt_CENT = - K20*CENT - K23*CENT + K32*PERI;
dxdt_PERI = K23*CENT - K32*PERI;

$OMEGA @annotated
ECL : 0.0484 : ETA on CL
EVC : 0.1110 : ETA on VC
EVP : 0.0451 : ETA on VP

$SIGMA @annotated
PROP: 0.05 : Proportional residual error

$TABLE
double no_obs_AUC = AMT/CLCOV;
double CP_no_RUV = CENT/VCCOV;
double CP = CP_no_RUV * (1 + PROP);
double CT = (PERI/VPCOV);

//prevent simulation of negative concentrations
int i = 0;
while(CP <0 && i < 100){
    simeps();
    CP = CP_no_RUV * (1 + PROP);
    ++i;
}

$capture
CP_no_RUV CP CT no_obs_AUC
