SharpClaw library files:

Core routines:

ClawData.f95   - Definition of clawdat and griddat data types
GlobalN.f95    - Sets shape of q, qold and qrk
qallocN.f95    - Allocate q, qold, and qrk
main.f95       - Outer routine to call sharpclaw and outN
sharpclaw.f95  - Inner loop between outputs
step.f95       - Take a single timestep (call, bcN, src, and stepN)
stepN.f95      - Homogeneous step for cons. law part
bcN.f95         - Boundary conditions (does this need to be dimensional?)
outN.f95       - Output (does this need to be dimensional?)
outn_hdf.f95   - HDF5 output

Reconstruction routines:

tvd2.f95    - TVD, 2nd order limiters
weno5.f95   - Fifth order WENO limiters
polyrecon.f95 - Non-limited reconstruction


Dummy routines to be overridden by user if needed:

b4step.f95
evec.f95
setaux.f95
setprob.f95
src.f95
tfluct.f95
