# GEOKEYFEM_HM in 2D
The numerical simulation code of [Qi ZHANG](https://qizhang94.github.io/). This work was supported by the **RGC Postdoctoral Fellowship Scheme** (RGC Ref. No. PDFS2223-5S04) and the **Start-up Fund for RAPs under the Strategic Hiring Scheme** (Grant No. P0043879).

[Buckley–Leverett Displacement](https://gitee.com/qzhang94/fem-for-buckley-leverett-displacement.git)

## Irregular updates
[A benchmark example](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/MAIN_Examples/main_2d_prandtl.m) that concerns the bearing capacity of foundation soil was added on 02/08/2023 (mm/dd/yyyy). The NS-FEM result could match perfectly with the Prandtl solution. The surface of discontinuity was also qualitatively correct. **Change the UMAT file (to D-P) in the assemble function**.

[Another example](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/MAIN_Examples/main_2d_Sloan_case.m) that tries to reproduce one case from the following paper was added on 06/03/2023 (mm/dd/yyyy). In order to match the reference result, the stabilization (tunning) parameter should be altered (a value of 1 would lead to inconsistent result). However, the **effective stress** distribution has some spurious oscillations, which is the weakness of SNS-FEM compared to standard FEM. By default, Mohr-Coulomb law should be used, despite that we use the Drucker-Prager model here.

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
The deformation equation is discretized by using the [Smoothed Finite Element Method](https://www.taylorfrancis.com/books/mono/10.1201/EBK1439820278/smoothed-finite-element-methods-liu-nguyen-trung). The direct nodal integration on the smoothing domain is also modified to **stabilized conforming nodal integration (SCNI)**. This is reflected in this [file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/assemble_stab.m).

The code is designed for hydromechanical coupling analysis (poromechanics). Therefore, the [pore pressure stabilization technique](https://doi.org/10.1016/j.cma.2008.05.015) is included. In addition, for simplicity, the biot coefficient is set to be `1` and Biot modulus is set to be `infinite`. It should be not too difficult to modify our code for a more general case by using relevant pre-defined matrices such as the `Mass matrix` "integral(N^T\*N)" from this [file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/pre_assemble_MassN.m).

You can run [this file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/main_2d_gujie.m) directly. The (max) strip load is 0.05 MPa applied in 100 seconds. The top surface is fully drained. In the main file, the most time-consuming part is the assembly process, whose function is given in this [file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/assemble_system.m). The [constitutive file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/DP_UMAT.m) is also called in the assembly process.

This [file](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/pre_assemble_BigK.m) *pre-aseembles* a **template** for the 2 by 2 block stiffness matrix. It makes the actual assemble [code](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/assemble_system.m) more concise. However, it only works for constant biot coefficient/tensor and mobility, i.e., some material properties will not change with (x,y,z).

The [assign_tractionBC2](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/assign_tractionBC2.m) is designed to calculate the equivalent nodal force vector "integral_line(N^T\*traction)" in FEM. The user only need to define the `traction_f` as a vectorial function of (x,y,z,time) (P.S. the element-wise multiplication `.*` should be adopted).

PFEM (Particle-FEM) feature has not been incorporated yet since it will affect the code structure.

The calculations of the equivalent plastic strain is based on the deviatoric plastic strain tensor.

## Output

- The undeformed mesh
- The deformed mesh
- The contour of horizontal and vertical displacements
- The contour of (effective) stress field (4 components: xx, yy, zz, xy)
- The contour of equivalent deviatoric plastic strain

