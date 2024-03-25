# About iMC code

iMC is in-house Monte Carlo neutron transport code, currently under development in RP&T lab in KAIST. 
iMC is based on Fortran90, supporting MPI- and OpenMP-based parallelization.

**Inhyung Kim** and **Hyeontae Kim** are initial contributors of the code.
**Inyup Kim**, **Taesuk Oh**, **Sunjoo Yoon**, and **Sangjin Lee** are current (2024.03.25) developers of the code.

The code supports following features
- Continuous and Multigroup transport
- Universe-based geometry building
- Steady-state MC
- OTF Doppler broadening
- MC-based depletion
- Transient with PCQS and DMC
- Finite Difference methods (pCMFD, pFMFD, and iDTMC)
- Unstructured mesh-based thermomechanics analysis
- ...
