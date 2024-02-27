# GEOKEYFEM_HM
The numerical simulation code of [Qi ZHANG](https://qizhang94.github.io/){:target="_blank"}. This work was supported by the **RGC Postdoctoral Fellowship Scheme** (RGC Ref. No. PDFS2223-5S04) and the **Start-up Fund for RAPs under the Strategic Hiring Scheme** (Grant No. P0043879).

**ALERT! ALERT! ALERT!** Please **DON'T** used the M-C UMAT code for **3D**, it will **NEVER WORK**! There are still many bugs! If you have to try 3D, modify the D-P UMAT code!

[Buckley–Leverett Displacement](https://gitee.com/qzhang94/fem-for-buckley-leverett-displacement.git){:target="_blank"}

## Irregular updates
[A benchmark example](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/main_2d_prandtl.m){:target="_blank"} that concerns the bearing capacity of foundation soil was added on 02/08/2023 (mm/dd/yyyy). The NS-FEM result could match perfectly with the Prandtl solution. The surface of discontinuity was also qualitatively correct. **Change the UMAT file (to D-P) in the assemble function**.

[Another example](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/main_2d_Sloan_case.m){:target="_blank"} that tries to reproduce one case from the following paper was added on 06/03/2023 (mm/dd/yyyy). In order to match the reference result, the stabilization (tunning) parameter should be altered (a value of 1 would lead to inconsistent result). However, the **effective stress** distribution has some spurious oscillations, which is the weakness of SNS-FEM compared to standard FEM. By default, Mohr-Coulomb law is used.

*Sloan, S.W. and Abbo, A.J. (1999), Biot consolidation analysis with automatic time stepping and error control Part 2: applications. Int. J. Numer. Anal. Meth. Geomech., 23: 493-529.*

```diff
! A new "GasHM" folder has been created,
! which describes how to incorporate gas into hydro-mechanical coupling FEM.
+ The folder contains an assembly function, a constitutive model for transversely isotropic geomaterial,
+ one verification example, and one application example.
@@ To run the example, you need to move files to the main directory. @@
```

### Reflection
In a nutshell, trade-off among **computational efficiency** (Highest: T6 FEM), **spurious oscillations** (numerical instability: NS-FEM based on T3 linear interpolation), and **overly rigid solution** (T3 FEM). SNS-FEM is in between.

## Functionality
The code simulates a contact problem between a rigid rectangular block with a Mohr-Coulomb soil by using the penalty method (**small deformation**). The deformation equation is discretized by using the [Smoothed Finite Element Method](https://www.taylorfrancis.com/books/mono/10.1201/EBK1439820278/smoothed-finite-element-methods-liu-nguyen-trung). The direct nodal integration on the smoothing domain is also modified to **stabilized conforming nodal integration (SCNI)**. This is reflected in this [file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/assemble_stab.m).


Although the current example only considers the deformation field, the code is designed for hydromechanical coupling analysis (poromechanics). Therefore, the [pore pressure stabilization technique](https://doi.org/10.1016/j.cma.2008.05.015) is included. In addition, for simplicity, the biot coefficient is set to be `1` and Biot modulus is set to be `infinite`. It should be not too difficult to modify our code for a more general case by using relevant pre-defined matrices such as the `Mass matrix` "integral(N^T\*N)" from this [file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/pre_assemble_MassN.m).


This is the [main file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/main_rigid_contact_prob.m). You can run it directly. In the main file, the most time-consuming part is the assembly process, whose function is given in this [file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/assemble_system.m). The [constitutive file](./MohrCoulomb_UMAT.m) is also called in the assembly process, as shown in the following code block:
```
[stress_new(:, ino), SDV_new(:, ino), cto] = MohrCoulomb_UMAT(0, Props(ino,:), stress(:, ino), strain_ino_new - strain_ino_old, SDV(:, ino));
```

This [file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/pre_assemble_BigK.m) *pre-aseembles* a **template** for the 2 by 2 block stiffness matrix. It makes the actual assemble [code](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/assemble_system.m) more concise. However, it only works for constant biot coefficient/tensor and mobility, i.e., some material properties will not change with (x,y,z).


The [assign_tractionBC2](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/assign_tractionBC2.m) is not used in this contact problem, while it is designed to calculate the equivalent nodal force vector "integral_line(N^T\*traction)" in FEM. The user only need to define the `traction_f` as a function of (x,y,z,time) (P.S. the element-wise multiplication `.*` should be adopted).

PFEM (Particle-FEM) feature has not been incorporated yet since it will affect the code structure.

The new main file considers a simple 2D Mohr-Coulomb soil consolidation problem. The geometry is the same as the contact problem. The (max) strip load is 0.05 MPa applied in 100 seconds. The top surface is fully drained.


## Advice

The code cannot be perfect. It is highly likely that non-convergence will happen if you **change some parameters**, especially when you want to push the rigid block forward. One scenario that would **partially work** is WHEN there is no friction, the initial location of the rigid block is {(0.29, 1) m, (0.51, 1) m}, the total displacement (10 time steps) is (0.3, -0.1) m, mesh size is 0.025 m, Mohr-Coulomb model is adopted, and `epsp` is 1. The PEEQ profile **cannot** match the standard FEM result.


However, if you consider friction such as making CFRI = 0.5, the original main code will not converge from the first time step. In that case, we have "by accidently" found a **new** scheme for updating stiffness matrix. This is given in the second main file with name `lucky`. We CANNOT guarantee that it will work for other examples.


The calculations of the equivalent plastic strain for D-P and M-C models are based on the deviatoric strain: **`depsp_d`**.


## Output
If you type `run main_rigid_contact_prob.m` in the MATLAB command window by using default parameters, you will get the [following sample output](./images/sample_out1_SNSFEM.txt){:target="_blank"} (note the *residual norm does not converge to zero, which is still under investigation*). Five figures will also be generated:
- The undeformed mesh
- The deformed mesh
- The contour of horizontal and vertical displacements
- The contour of (effective) stress field (4 components: xx, yy, zz, xy)
- The contour of equivalent deviatoric plastic strain

An overview of the graphical result (**NOT reference result!**) could be found in the following [Tencent document](https://docs.qq.com/doc/DZUlGcndWaVBud0NP).


[Sample output](./images/sample_out2_SNSFEM.txt){:target="_blank"} from MATLAB command window of the Sloan and Abbo (1999) case (Mohr-Coulomb plasticity).


## To do list
- Divide the contact surface into many small segments to see the effectiveness.
- Consider objective stress rate (by spin tensor **w**).
- Re-meshing strategy for large deformation simulation.

## Important computer folders
By default, computer folders are based on the LG laptop.

#### Key
- C:\Users\zq112\Desktop\ (COMSOL notes)
- E:\Numerical simulation\COMSOL_NEWsimulations\Zhang2021_CMAT_2nd_example.mph (Published paper that uses COMSOL)
- E:\Postdoc_PolyU\Coder_MEX\GEOKEYFEM_HM (Github repo)
- C:\Users\zq112\OneDrive\REMOTE_SYNC\F\CSE583_Analytical and Numerical Methods in Geotechnical Engineering (Useful functions such as the *anisotropic elasticity* or *EVP* or *isotropic function*)
- Overleaf CSE583 Notes about MCC

#### On-going paper and revision
- **C:\Users\CEE\Desktop**\kapp\main_test_validate_coal_Aug08.m + main_test_for_COMSOL.m (Remote PolyU computer)
- **C:\Users\CEE\Desktop**\kapp\New_simulations\Gas_transfer_PLA_ksdiff(surf_dif)_perm_bedding_BN=BP_add_newadsorption_PLASTIC.mph (Remote PolyU computer) (Permeability is isotropic, but elastic tensor isn't, theta = pi/6)
- √ C:\Users\zq112\OneDrive\期刊论文修改发表2022\JMPS_or_IJNAG\Revision (Compare adsorbed mass functions)
- √ C:\GEOKEYFEM_HM\gas_data (ASUS Gaming Laptop)


================

- √ C:\Users\zq112\OneDrive\REMOTE_SYNC\E\Postdoc_PolyU\Shale_gas_on_2022_not_published\Verification\New_case_anisotropy_3D.mph (For verification)
- √ C:\Users\zq112\OneDrive\REMOTE_SYNC\E\Postdoc_PolyU\Shale_gas_on_2022_not_published\POROELASTIC_CONSTANTS.xlsx
- √ C:\Users\zq112\OneDrive\REMOTE_SYNC\E\Postdoc_PolyU\Shale_gas_on_2022_not_published\Revised_results.ai
- **C:\Users\CEE\Desktop**\3D_reservoir_COMSOL (Remote PolyU computer, see Excel subfolder)
- **C:\Users\CEE\Desktop** Verification_3D1_ANISOTROPIC_horizontal_!!!!!!!!!!!!!!!!!!!!!_IMPORTANT.mph on the Desktop folder of the remote PolyU computer

#### Finite strain solid mechanics
- C:\Users\zq112\OneDrive\REMOTE_SYNC\Nonlinear FEA & C. Linder
- C:\Users\zq112\OneDrive\REMOTE_SYNC\E\Postdoc_PolyU\NFEA
- C:\Users\zq112\OneDrive\REMOTE_SYNC\E\Postdoc_PolyU\NFEA\T3_Solid_To_ZQ

#### Group
- C:\Users\zq112\Desktop\SNS-FEM-MCC-NEW3 (MCC + large deformation remesh)
- C:\Users\zq112\OneDrive\REMOTE_SYNC\E\Postdoc_PolyU\Zeyu_WANG\FEMCON_T3T36 (Standard FEM, coded by Ze-Yu WANG)
- C:\Users\zq112\OneDrive\REMOTE_SYNC\E\FEM_HM_T6T3 (Standard FEM, coded by Xianhan WU)

#### Hydrate
- C:\Users\zq112\Desktop\1D_hydrate\reproduce_Klar.m
- C:\Users\zq112\Desktop\1D_hydrate\plot_compare.m
- C:\Users\CEE\Desktop\1D_hydrate\main_solve_thermo_Masuda2_true_sequential.m (Remote PolyU computer)
- C:\Users\zq112\OneDrive\REMOTE_SYNC\F\CSE583_Analytical and Numerical Methods in Geotechnical Engineering\DP_drained_or_undrained (optional)\CDtest_Hydrate_Sh.m
