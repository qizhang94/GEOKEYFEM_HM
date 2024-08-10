## Advice

The code cannot be perfect. It is highly likely that non-convergence will happen if you **change some parameters**, especially when you want to push the rigid block forward. One scenario that would **partially work** is WHEN there is no friction, the initial location of the rigid block is {(0.29, 1) m, (0.51, 1) m}, the total displacement (10 time steps) is (0.3, -0.1) m, mesh size is 0.025 m, Mohr-Coulomb model is adopted, and `epsp` is 1. The PEEQ profile **cannot** match the standard FEM result.


However, if you consider friction such as making CFRI = 0.5, the original main code will not converge from the first time step. In that case, we have "by accidently" found a **new** scheme for updating stiffness matrix. This is given in the second main file with name `lucky`. We CANNOT guarantee that it will work for other examples.


 for D-P and M-C models are based on the deviatoric strain: **`depsp_d`**. Note in these models, the plastic strain tensor is ususally **NOT** stored using the Voigt notation.


## Output
If you type `run main_rigid_contact_prob.m` in the MATLAB command window by using default parameters, you will get the [following sample output](./images/sample_out1_SNSFEM.txt) (note the *residual norm does not converge to zero, which is still under investigation*). Five figures will also be generated:
- The undeformed mesh
- The deformed mesh
- The contour of horizontal and vertical displacements
- The contour of (effective) stress field (4 components: xx, yy, zz, xy)
- The contour of equivalent deviatoric plastic strain

An overview of the graphical result (**NOT reference result!**) could be found in the following [Tencent document](https://docs.qq.com/doc/DZUlGcndWaVBud0NP).


[Sample output](./images/sample_out2_SNSFEM.txt) from MATLAB command window of the Sloan and Abbo (1999) case (Mohr-Coulomb plasticity).


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