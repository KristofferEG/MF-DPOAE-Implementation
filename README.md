# MF-DPOAE-Implementation
This is a repository containing code for a Multi Frequency DPOAE implementation for the Interacoustics Titan. 

The files of note are the following:
**- measMFDPOA.m** \\
  The main MF–DPOAE code.
**-  MFDPOAEtestBattery.mlx**
  The code used to automate and visualize testing. It also contains the different DPOAE/MF-DPOAE variations at the top, with the current settings in the code.
**- generateStimulus.m**
  The code which generates the stimulus/primary frequencies ( f1, f2).
**- dp_frequency_test2.mlx**
  The full frequency pairing algorithm. Protects from overlap with most DPs.
**- dp_frequency_test2Lax.mlx**
  The Relaxed frequency pairing algorithm. Protects from overlap with fewer DPs. 

The rest of the files are related to the report work performed as part of its development.

The different variations of DPOAE/MF-DPOAE which can be performed using this is:
- 6 frequency standard DPOAE
- 12 frequency standard DPOAE
- 6 frequency manual minimum octave spacing
- 12 frequency manual minimum octave spacing
- 6 frequency algorithm based spacing
- 12 frequency algorithm based spacing
- 6 frequency relaxed algorithm based spacing
- 12 frequency relaxed algorithm based spacing
