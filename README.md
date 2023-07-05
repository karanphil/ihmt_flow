# ihMT pipeline
===================

Flow to process ihMT images. This flow generated four myelin indices maps from the inhomogeneous Magnetization Transfer (ihMT) images. 

Input Data recommendation:
  - it is recommended to use dcm2niix (v1.0.20200331) to convert data
    https://github.com/rordenlab/dcm2niix/releases/tag/v1.0.20200331
  - dcm2niix conversion will create all echo files for each contrast and
    corresponding json files
  - all input must have a matching json file with the same filename
  - all contrasts must have a same number of echoes and coregistered
    between them before running the script.
  - Mask must be coregistered to the echo images
  - ANTs can be used for the registration steps (http://stnava.github.io/ANTs/)


A full description is available in this reference: https://onlinelibrary.wiley.com/doi/10.1002/hbm.26310. 
The supplementary data provide the ihMT maps computations in more detail. 



Requirements
------------

- [Nextflow](https://www.nextflow.io)
- [scilpy](https://github.com/scilus/scilpy)
- [ants](https://github.com/ANTsX/ANTs)



Singularity
-----------

If you are on Linux, we recommend using the Singularity container to run ihMT_flow

Prebuild Singularity images are available here:

[https://scil.usherbrooke.ca/pages/containers/](https://scil.usherbrooke.ca/pages/containers/)

FOR DEVELOPERS: The containers repository is available here:
[containers-scilus](https://github.com/scilus/containers-scilus)



Steps
-----

- Coregistration of images (ANTs)
- Brain extraction (FSL bet)
- Segmentation (ANTs)
- Compute ihMT maps (Scilpy)
- Registration to diffusion space using Register T1 of Tractoflow (ANTs)



Output
------

- ihMT dR1 saturation (ihMTsat)
- ihMT ratio (ihMTR)
- MT saturation (MTsat)
- MT Ratio (MTR)

Myelin-sensitive maps in native space: Contrasts_ihMT_maps and ihMT_native_maps folders.

Myelin-sensitive maps in diffusion space: Register_contrast_maps and Register_ihMT_maps folders.

** Update for PK (main_pk_analysis.nf): 
Include Contrast maps in Diffusion space as output.

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`


References
----------

```
Varma G, Girard OM, Prevost VH, Grant AK, Duhamel G, Alsop DC.
Interpretation of magnetization transfer from inhomogeneously broadened lines
(ihMT) in tissues as a dipolar order effect within motion restricted molecules.
Journal of Magnetic Resonance. 1 nov 2015;260:67-76.

Manning AP, Chang KL, MacKay AL, Michal CA. The physical mechanism of
"inhomogeneous" magnetization transfer MRI. Journal of Magnetic Resonance.
1 janv 2017;274:125-36.

Helms G, Dathe H, Kallenberg K, Dechent P. High-resolution maps of
magnetization transfer with inherent correction for RF inhomogeneity
and T1 relaxation obtained from 3D FLASH MRI. Magnetic Resonance in Medicine.
2008;60(6):1396-407.
"""
```
