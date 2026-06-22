# purkinje-ventricular-junction
computational simulation of EAD at the purkinje ventricular junction

Here is the plain text version:
Purkinje-Ventricular Junction Simulation
This directory contains Fortran source code for simulating calcium-mediated retrograde propagation at a simplified Purkinje-ventricular junction, or PVJ.
The model couples stochastic intracellular calcium cycling to membrane electrophysiology in two connected tissue regions:
A narrow Purkinje-like fiber, represented as a 100 x 6 grid.
A ventricular tissue sheet, represented as a 100 x 52 grid in the current baseline code.
The code is intended for research simulations and figure generation for the Purkinje propagation project. Before changing parameters or model equations, users should read the accompanying manuscript and online supplement:
purkinje_propagation.pdf
online-supplement.pdf
Source Files
space5.f
Main simulation program. Defines the geometry, pacing protocol, model parameters, main time loop, tissue updates, and output files.
jun3.f
Initialization routines and voltage diffusion/coupling between the Purkinje-like strip and ventricular sheet.
avcurrents7.f
Membrane ionic current routines, including INa, IKr, IKs, IK1, IKur, Ito, and INaK.
acur6.f
Calcium cycling support routines, stochastic spark evolution, buffering, NCX, L-type Ca current, and the L-type Ca channel Markov model.
The program must be compiled from all four source files together.
Required Compiler
Use the Intel Fortran compiler ifx.
On Windows, open an Intel oneAPI command prompt or a terminal where ifx is available on the PATH.
Check that the compiler is available with:
ifx --version
Compile
From the code directory, compile with:
ifx space5.f jun3.f avcurrents7.f acur6.f /O2 /exe:space5.exe
For a debug build with runtime checks, use:
ifx space5.f jun3.f avcurrents7.f acur6.f /debug /check:all /traceback /exe:space5_debug.exe
The optimized build is appropriate for production simulations. The debug build is useful after code changes, especially if array bounds or uninitialized values are suspected.
Run
Run the executable from the same directory:
.\space5.exe
The program writes output files into the current working directory. To keep runs organized, it is best to make a separate folder for each run, copy the executable into that folder, and run it there.
Baseline Case
The current baseline geometry is:
Lx = 100
Ly = 6
Lx1 = 100
Ly1 = 52
This corresponds to a 100 x 6 Purkinje-like fiber connected to a 100 x 52 ventricular sheet.
Important baseline parameters include:
bcl = 640.0
dt = 0.1 / 1.50
gica = 1.5
gicaz = 1.8
pbx = 0.3
pbxz = 0.3
Here gica is the L-type Ca conductance scaling in the Purkinje-like strip, while gicaz is the corresponding value in the ventricular sheet. The larger ventricular value makes the sheet EAD-prone.
Main Output Files
vs10.dat
Time and voltage at v(10,2) in the Purkinje-like strip.
vs100.dat
Time and voltage at v(100,2), near the junction end of the strip.
vs190.dat
Time and voltage at vz(80,2) in the ventricular sheet.
icas10.dat
Time and a Ca-current diagnostic at wica(20,2).
Space-Time Output Files
The code also writes numbered files along the centerline of the combined geometry. The first Lx points correspond to the Purkinje-like strip, and the next Lx1 points correspond to the ventricular sheet.
vb*.dat contains voltage along the combined centerline.
ik*.dat contains a potassium-current diagnostic.
icaz*.dat contains a Ca-current diagnostic.
inaz*.dat contains a sodium-current diagnostic.
po*.dat contains a sodium activation proxy.
poinf*.dat contains a sodium availability proxy.
pox*.dat contains a product-like sodium recruitment proxy.
Each file contains two columns:
position_index   value
These files are used to build space-time plots like those in the manuscript.
Reproducibility Notes
This code includes stochastic calcium spark evolution through calls to random_number. For exact reproducibility across runs, a fixed random seed should be added and recorded. Until then, repeated runs may not be bit-for-bit identical.
For every simulation used in a paper figure, record:
date of run
code version or Git commit
compiler and compiler version
values of Lx, Ly, Lx1, and Ly1
bcl and dt
gica and gicaz
pbx and pbxz
any other changed parameters
output folder
plotting script used
brief note on whether retrograde propagation succeeded or failed
Important Cautions
Small parameter changes can shift the transition between propagation failure and successful retrograde excitation.
The current code uses Ly1 = 52, a near-threshold successful case in the manuscript.
The supplement should be checked against the code before publication figures are finalized. For example, the code currently uses nit = 6000.
The routine markovz exists in acur6.f, but space5.f currently calls markov for both tissue regions.
The diffusion/coupling routine Euler_Forward is called twice in each main time step. Do not change this unless deliberately testing the numerical method.
Quick Start
Compile:
ifx space5.f jun3.f avcurrents7.f acur6.f /O2 /exe:space5.exe
Run:
.\space5.exe
After the run completes, inspect vs10.dat, vs100.dat, and vs190.dat. Then use the numbered vb*.dat, icaz*.dat, and inaz*.dat files to generate space-time plots.
