# GEOKEYFEM_HM
The numerical simulation code of [Qi ZHANG](https://qizhang94.github.io/). This work was supported by the **RGC Postdoctoral Fellowship Scheme** (RGC Ref. No. PDFS2223-5S04) and the **Start-up Fund for RAPs under the Strategic Hiring Scheme** (Grant No. P0043879).

**ALERT! ALERT! ALERT!** Please **DON'T** used the M-C UMAT code for **3D**, it will **NEVER WORK**! There are still many bugs! If you have to try 3D, modify the D-P UMAT code!

## Irregular updates
[A benchmark example](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/main_2d_prandtl.m) that concerns the bearing capacity of foundation soil was added on 02/08/2023 (mm/dd/yyyy). The NS-FEM result could match perfectly with the Prandtl solution. The surface of discontinuity was also qualitatively correct. **Change the UMAT file (to D-P) in the assemble function**.

[Another example](https://github.com/qizhang94/GEOKEYFEM_HM/blob/main/main_2d_Sloan_case.m) that tries to reproduce one case from the following paper was added on 06/03/2023 (mm/dd/yyyy). In order to match the reference result, the stabilization (tunning) parameter should be altered (a value of 1 would lead to inconsistent result). However, the **effective stress** distribution has some spurious oscillations, which is the weakness of SNS-FEM compared to standard FEM. By default, Mohr-Coulomb law is used.

*Sloan, S.W. and Abbo, A.J. (1999), Biot consolidation analysis with automatic time stepping and error control Part 2: applications. Int. J. Numer. Anal. Meth. Geomech., 23: 493-529.*

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


However, if you consider friction such as making CFRI = 0.5, the original main code will not converge from the first time step. In that case, we have found a **new** scheme for updating stiffness matrix "by accident". This is given in the second main file with name `lucky`. We CANNOT guarantee that it will work for other examples.


The calculations of the equivalent plastic strain for D-P and M-C models are based on the deviatoric strain: **`depsp_d`**.


## Output
If you type `run main_rigid_contact_prob.m` in the MATLAB command window by using default parameters, you will get the following sample output (note the *residual norm does not converge to zero, which is still under investigation*). Five figures will also be generated:
- The undeformed mesh
- The deformed mesh
- The contour of horizontal and vertical displacements
- The contour of (effective) stress field (4 components: xx, yy, zz, xy)
- The contour of equivalent deviatoric plastic strain

An overview of the graphical result could be found in the following [Tencent document](https://docs.qq.com/doc/DZUlGcndWaVBud0NP).

```
LOAD STEP = 1; TIME = 8.64:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 1.2921E-10; RES = 5.6370E-05
LOAD STEP = 2; TIME = 17.28:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 2.4044E-10; RES = 1.0359E-04
LOAD STEP = 3; TIME = 25.92:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 1.7993E-04; RES = 1.4820E-04
LOAD STEP = 4; TIME = 34.56:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 3.7912E-03; RES = 1.8788E-04
	 NEWTON = 3; ERROR = 1.3685E-04; RES = 1.8766E-04
LOAD STEP = 5; TIME = 43.20:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 9.3549E-03; RES = 2.2662E-04
	 NEWTON = 3; ERROR = 7.9673E-04; RES = 2.2310E-04
LOAD STEP = 6; TIME = 51.84:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 2.3144E-02; RES = 2.5181E-04
	 NEWTON = 3; ERROR = 3.9380E-03; RES = 2.4277E-04
	 NEWTON = 4; ERROR = 8.6805E-04; RES = 2.4141E-04
LOAD STEP = 7; TIME = 60.48:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 4.0250E-02; RES = 2.6959E-04
	 NEWTON = 3; ERROR = 1.0140E-02; RES = 2.5473E-04
	 NEWTON = 4; ERROR = 2.8344E-03; RES = 2.5249E-04
	 NEWTON = 5; ERROR = 7.6317E-04; RES = 2.5204E-04
LOAD STEP = 8; TIME = 69.12:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 5.6692E-02; RES = 2.7744E-04
	 NEWTON = 3; ERROR = 2.0020E-02; RES = 2.5996E-04
	 NEWTON = 4; ERROR = 7.0351E-03; RES = 2.5681E-04
	 NEWTON = 5; ERROR = 2.4637E-03; RES = 2.5616E-04
	 NEWTON = 6; ERROR = 9.2800E-04; RES = 2.5595E-04
LOAD STEP = 9; TIME = 77.76:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 1.2716E-01; RES = 2.9601E-04
	 NEWTON = 3; ERROR = 4.6294E-02; RES = 2.6010E-04
	 NEWTON = 4; ERROR = 1.9707E-02; RES = 2.5343E-04
	 NEWTON = 5; ERROR = 7.3479E-03; RES = 2.5172E-04
	 NEWTON = 6; ERROR = 3.2769E-03; RES = 2.5123E-04
	 NEWTON = 7; ERROR = 1.3195E-03; RES = 2.5108E-04
	 NEWTON = 8; ERROR = 6.0660E-04; RES = 2.5103E-04
LOAD STEP = 10; TIME = 86.40:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0000E+06
	 NEWTON = 2; ERROR = 2.0702E-01; RES = 2.8667E-04
	 NEWTON = 3; ERROR = 8.3940E-02; RES = 2.5344E-04
	 NEWTON = 4; ERROR = 3.6587E-02; RES = 2.4616E-04
	 NEWTON = 5; ERROR = 1.5123E-02; RES = 2.4407E-04
	 NEWTON = 6; ERROR = 6.9927E-03; RES = 2.4348E-04
	 NEWTON = 7; ERROR = 3.1249E-03; RES = 2.4328E-04
	 NEWTON = 8; ERROR = 1.5108E-03; RES = 2.4321E-04
	 NEWTON = 9; ERROR = 7.1468E-04; RES = 2.4318E-04
```

Sample output from MATLAB command window of the Sloan and Abbo (1999) case.

```
LOAD STEP = 1; TIME = 54912.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.1227E-03
	 NEWTON = 2; ERROR = 8.1434E-16; RES = 2.2457E-18
LOAD STEP = 2; TIME = 109824.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.1535E-03
	 NEWTON = 2; ERROR = 9.9107E-16; RES = 3.2208E-18
LOAD STEP = 3; TIME = 164736.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.1675E-03
	 NEWTON = 2; ERROR = 5.9210E-16; RES = 3.1756E-18
LOAD STEP = 4; TIME = 219648.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.1770E-03
	 NEWTON = 2; ERROR = 1.3337E-15; RES = 4.2471E-18
LOAD STEP = 5; TIME = 274560.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.1843E-03
	 NEWTON = 2; ERROR = 8.5926E-16; RES = 5.5007E-18
LOAD STEP = 6; TIME = 329472.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.1900E-03
	 NEWTON = 2; ERROR = 9.9942E-16; RES = 5.1807E-18
LOAD STEP = 7; TIME = 384384.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.1947E-03
	 NEWTON = 2; ERROR = 6.9998E-16; RES = 6.6937E-18
LOAD STEP = 8; TIME = 439296.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.1987E-03
	 NEWTON = 2; ERROR = 1.1167E-15; RES = 6.2070E-18
LOAD STEP = 9; TIME = 494208.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2021E-03
	 NEWTON = 2; ERROR = 1.0974E-15; RES = 8.1017E-18
LOAD STEP = 10; TIME = 549120.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2050E-03
	 NEWTON = 2; ERROR = 1.4282E-15; RES = 1.2130E-17
LOAD STEP = 11; TIME = 604032.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2075E-03
	 NEWTON = 2; ERROR = 1.6619E-15; RES = 1.1986E-17
LOAD STEP = 12; TIME = 658944.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2097E-03
	 NEWTON = 2; ERROR = 1.6714E-15; RES = 1.2456E-17
LOAD STEP = 13; TIME = 713856.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2116E-03
	 NEWTON = 2; ERROR = 1.9280E-15; RES = 1.4375E-17
LOAD STEP = 14; TIME = 768768.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2134E-03
	 NEWTON = 2; ERROR = 2.3475E-15; RES = 1.3598E-17
LOAD STEP = 15; TIME = 823680.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2149E-03
	 NEWTON = 2; ERROR = 1.8821E-15; RES = 1.6759E-17
LOAD STEP = 16; TIME = 878592.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2163E-03
	 NEWTON = 2; ERROR = 5.9282E-02; RES = 1.4077E-04
	 NEWTON = 3; ERROR = 1.7251E-02; RES = 5.6235E-05
	 NEWTON = 4; ERROR = 5.3367E-03; RES = 1.6574E-05
	 NEWTON = 5; ERROR = 1.7135E-03; RES = 5.4498E-06
	 NEWTON = 6; ERROR = 5.5869E-04; RES = 1.7930E-06
LOAD STEP = 17; TIME = 933504.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2248E-03
	 NEWTON = 2; ERROR = 1.8878E-01; RES = 4.4363E-04
	 NEWTON = 3; ERROR = 8.7506E-02; RES = 2.2014E-04
	 NEWTON = 4; ERROR = 4.0936E-02; RES = 1.1133E-04
	 NEWTON = 5; ERROR = 1.9256E-02; RES = 5.1410E-05
	 NEWTON = 6; ERROR = 9.1063E-03; RES = 2.4224E-05
	 NEWTON = 7; ERROR = 4.3156E-03; RES = 1.1464E-05
	 NEWTON = 8; ERROR = 2.0473E-03; RES = 5.4355E-06
	 NEWTON = 9; ERROR = 9.7172E-04; RES = 2.5791E-06
LOAD STEP = 18; TIME = 988416.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.2678E-03
	 NEWTON = 2; ERROR = 3.0228E-01; RES = 5.6179E-04
	 NEWTON = 3; ERROR = 1.3697E-01; RES = 3.9522E-04
	 NEWTON = 4; ERROR = 6.7840E-02; RES = 1.9324E-04
	 NEWTON = 5; ERROR = 3.4228E-02; RES = 1.0421E-04
	 NEWTON = 6; ERROR = 1.7583E-02; RES = 5.3954E-05
	 NEWTON = 7; ERROR = 9.0712E-03; RES = 2.8239E-05
	 NEWTON = 8; ERROR = 4.6987E-03; RES = 1.4671E-05
	 NEWTON = 9; ERROR = 2.4369E-03; RES = 7.6336E-06
	 NEWTON = 10; ERROR = 1.2651E-03; RES = 3.9665E-06
	 NEWTON = 11; ERROR = 6.5700E-04; RES = 2.0615E-06
LOAD STEP = 19; TIME = 1043328.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.3448E-03
	 NEWTON = 2; ERROR = 3.2826E-01; RES = 6.6891E-04
	 NEWTON = 3; ERROR = 1.5563E-01; RES = 4.1836E-04
	 NEWTON = 4; ERROR = 7.8011E-02; RES = 2.1222E-04
	 NEWTON = 5; ERROR = 4.0110E-02; RES = 1.1683E-04
	 NEWTON = 6; ERROR = 2.1013E-02; RES = 6.2496E-05
	 NEWTON = 7; ERROR = 1.1083E-02; RES = 3.3580E-05
	 NEWTON = 8; ERROR = 5.8728E-03; RES = 1.7932E-05
	 NEWTON = 9; ERROR = 3.1184E-03; RES = 9.5720E-06
	 NEWTON = 10; ERROR = 1.6579E-03; RES = 5.1018E-06
	 NEWTON = 11; ERROR = 8.8195E-04; RES = 2.7183E-06
LOAD STEP = 20; TIME = 1098240.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.4358E-03
	 NEWTON = 2; ERROR = 3.8215E-01; RES = 7.1713E-04
	 NEWTON = 3; ERROR = 1.8337E-01; RES = 5.1262E-04
	 NEWTON = 4; ERROR = 9.9352E-02; RES = 2.6585E-04
	 NEWTON = 5; ERROR = 5.5052E-02; RES = 1.6244E-04
	 NEWTON = 6; ERROR = 3.1300E-02; RES = 9.4258E-05
	 NEWTON = 7; ERROR = 1.7947E-02; RES = 5.5406E-05
	 NEWTON = 8; ERROR = 1.0358E-02; RES = 3.2275E-05
	 NEWTON = 9; ERROR = 5.9966E-03; RES = 1.8824E-05
	 NEWTON = 10; ERROR = 3.5623E-03; RES = 1.0949E-05
	 NEWTON = 11; ERROR = 2.1007E-03; RES = 6.5585E-06
	 NEWTON = 12; ERROR = 1.2466E-03; RES = 3.9087E-06
	 NEWTON = 13; ERROR = 7.4000E-04; RES = 2.3298E-06
LOAD STEP = 21; TIME = 1153152.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.5380E-03
	 NEWTON = 2; ERROR = 4.2060E-01; RES = 7.8419E-04
	 NEWTON = 3; ERROR = 2.0295E-01; RES = 6.1491E-04
	 NEWTON = 4; ERROR = 1.1213E-01; RES = 3.4064E-04
	 NEWTON = 5; ERROR = 6.3708E-02; RES = 2.1083E-04
	 NEWTON = 6; ERROR = 3.7118E-02; RES = 1.2587E-04
	 NEWTON = 7; ERROR = 2.1843E-02; RES = 7.5971E-05
	 NEWTON = 8; ERROR = 1.2948E-02; RES = 4.5563E-05
	 NEWTON = 9; ERROR = 7.7035E-03; RES = 2.7341E-05
	 NEWTON = 10; ERROR = 4.5943E-03; RES = 1.6382E-05
	 NEWTON = 11; ERROR = 2.7437E-03; RES = 9.8132E-06
	 NEWTON = 12; ERROR = 1.6399E-03; RES = 5.8757E-06
	 NEWTON = 13; ERROR = 9.8066E-04; RES = 3.5175E-06
LOAD STEP = 22; TIME = 1208064.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.7307E-03
	 NEWTON = 2; ERROR = 4.3036E-01; RES = 8.5023E-04
	 NEWTON = 3; ERROR = 2.2446E-01; RES = 6.2120E-04
	 NEWTON = 4; ERROR = 1.2549E-01; RES = 3.9631E-04
	 NEWTON = 5; ERROR = 7.3868E-02; RES = 2.5254E-04
	 NEWTON = 6; ERROR = 4.4513E-02; RES = 1.5829E-04
	 NEWTON = 7; ERROR = 2.7164E-02; RES = 9.9154E-05
	 NEWTON = 8; ERROR = 1.6705E-02; RES = 6.1924E-05
	 NEWTON = 9; ERROR = 1.0320E-02; RES = 3.8637E-05
	 NEWTON = 10; ERROR = 6.3937E-03; RES = 2.4083E-05
	 NEWTON = 11; ERROR = 3.9681E-03; RES = 1.5004E-05
	 NEWTON = 12; ERROR = 2.4654E-03; RES = 9.3446E-06
	 NEWTON = 13; ERROR = 1.5328E-03; RES = 5.8184E-06
	 NEWTON = 14; ERROR = 9.5333E-04; RES = 3.6224E-06
LOAD STEP = 23; TIME = 1262976.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 1.9366E-03
	 NEWTON = 2; ERROR = 4.7899E-01; RES = 9.3927E-04
	 NEWTON = 3; ERROR = 2.4280E-01; RES = 7.4683E-04
	 NEWTON = 4; ERROR = 1.3825E-01; RES = 4.5767E-04
	 NEWTON = 5; ERROR = 8.2084E-02; RES = 2.9585E-04
	 NEWTON = 6; ERROR = 4.9993E-02; RES = 1.8755E-04
	 NEWTON = 7; ERROR = 3.0864E-02; RES = 1.1926E-04
	 NEWTON = 8; ERROR = 1.9219E-02; RES = 7.5555E-05
	 NEWTON = 9; ERROR = 1.2437E-02; RES = 4.7717E-05
	 NEWTON = 10; ERROR = 7.8299E-03; RES = 3.0868E-05
	 NEWTON = 11; ERROR = 4.9830E-03; RES = 1.9612E-05
	 NEWTON = 12; ERROR = 3.1685E-03; RES = 1.2551E-05
	 NEWTON = 13; ERROR = 2.0180E-03; RES = 8.0116E-06
	 NEWTON = 14; ERROR = 1.2859E-03; RES = 5.1154E-06
	 NEWTON = 15; ERROR = 8.1972E-04; RES = 3.2648E-06
LOAD STEP = 24; TIME = 1317888.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 2.1283E-03
	 NEWTON = 2; ERROR = 5.2344E-01; RES = 1.0167E-03
	 NEWTON = 3; ERROR = 2.5222E-01; RES = 8.6282E-04
	 NEWTON = 4; ERROR = 1.4478E-01; RES = 5.0159E-04
	 NEWTON = 5; ERROR = 8.6161E-02; RES = 3.2952E-04
	 NEWTON = 6; ERROR = 5.2855E-02; RES = 2.0926E-04
	 NEWTON = 7; ERROR = 3.2865E-02; RES = 1.3424E-04
	 NEWTON = 8; ERROR = 2.0626E-02; RES = 8.5710E-05
	 NEWTON = 9; ERROR = 1.3014E-02; RES = 5.4742E-05
	 NEWTON = 10; ERROR = 8.2402E-03; RES = 3.4923E-05
	 NEWTON = 11; ERROR = 5.2286E-03; RES = 2.2271E-05
	 NEWTON = 12; ERROR = 3.3223E-03; RES = 1.4196E-05
	 NEWTON = 13; ERROR = 2.1128E-03; RES = 9.0471E-06
	 NEWTON = 14; ERROR = 1.3444E-03; RES = 5.7645E-06
	 NEWTON = 15; ERROR = 8.5576E-04; RES = 3.6725E-06
LOAD STEP = 25; TIME = 1372800.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 2.3095E-03
	 NEWTON = 2; ERROR = 5.2305E-01; RES = 1.1920E-03
	 NEWTON = 3; ERROR = 2.5959E-01; RES = 9.2011E-04
	 NEWTON = 4; ERROR = 1.5209E-01; RES = 5.4034E-04
	 NEWTON = 5; ERROR = 9.1945E-02; RES = 3.5685E-04
	 NEWTON = 6; ERROR = 5.7951E-02; RES = 2.3151E-04
	 NEWTON = 7; ERROR = 3.7969E-02; RES = 1.5905E-04
	 NEWTON = 8; ERROR = 2.4400E-02; RES = 1.0301E-04
	 NEWTON = 9; ERROR = 1.6125E-02; RES = 6.9443E-05
	 NEWTON = 10; ERROR = 1.0507E-02; RES = 4.5426E-05
	 NEWTON = 11; ERROR = 6.9587E-03; RES = 3.0377E-05
	 NEWTON = 12; ERROR = 4.5663E-03; RES = 1.9951E-05
	 NEWTON = 13; ERROR = 3.0250E-03; RES = 1.3283E-05
	 NEWTON = 14; ERROR = 1.9921E-03; RES = 8.7458E-06
	 NEWTON = 15; ERROR = 1.3193E-03; RES = 5.8074E-06
	 NEWTON = 16; ERROR = 8.7049E-04; RES = 3.8301E-06
LOAD STEP = 26; TIME = 1427712.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 2.5385E-03
	 NEWTON = 2; ERROR = 5.4827E-01; RES = 1.1745E-03
	 NEWTON = 3; ERROR = 2.7754E-01; RES = 1.0086E-03
	 NEWTON = 4; ERROR = 1.6409E-01; RES = 6.1977E-04
	 NEWTON = 5; ERROR = 1.0309E-01; RES = 4.3100E-04
	 NEWTON = 6; ERROR = 6.6641E-02; RES = 2.9251E-04
	 NEWTON = 7; ERROR = 4.4137E-02; RES = 2.0343E-04
	 NEWTON = 8; ERROR = 2.9435E-02; RES = 1.3970E-04
	 NEWTON = 9; ERROR = 1.9865E-02; RES = 9.6442E-05
	 NEWTON = 10; ERROR = 1.3430E-02; RES = 6.6107E-05
	 NEWTON = 11; ERROR = 9.1332E-03; RES = 4.5428E-05
	 NEWTON = 12; ERROR = 6.2132E-03; RES = 3.1098E-05
	 NEWTON = 13; ERROR = 4.2391E-03; RES = 2.1321E-05
	 NEWTON = 14; ERROR = 2.8920E-03; RES = 1.4587E-05
	 NEWTON = 15; ERROR = 1.9760E-03; RES = 9.9895E-06
	 NEWTON = 16; ERROR = 1.3499E-03; RES = 6.8331E-06
	 NEWTON = 17; ERROR = 9.2291E-04; RES = 4.6767E-06
LOAD STEP = 27; TIME = 1482624.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 2.7922E-03
	 NEWTON = 2; ERROR = 5.9830E-01; RES = 1.2990E-03
	 NEWTON = 3; ERROR = 2.9925E-01; RES = 1.1797E-03
	 NEWTON = 4; ERROR = 1.8145E-01; RES = 7.3203E-04
	 NEWTON = 5; ERROR = 1.1442E-01; RES = 5.2286E-04
	 NEWTON = 6; ERROR = 7.4730E-02; RES = 3.6161E-04
	 NEWTON = 7; ERROR = 4.9929E-02; RES = 2.5330E-04
	 NEWTON = 8; ERROR = 3.3652E-02; RES = 1.7549E-04
	 NEWTON = 9; ERROR = 2.2922E-02; RES = 1.2217E-04
	 NEWTON = 10; ERROR = 1.5662E-02; RES = 8.4577E-05
	 NEWTON = 11; ERROR = 1.0755E-02; RES = 5.8693E-05
	 NEWTON = 12; ERROR = 7.3950E-03; RES = 4.0608E-05
	 NEWTON = 13; ERROR = 5.0972E-03; RES = 2.8133E-05
	 NEWTON = 14; ERROR = 3.5151E-03; RES = 1.9458E-05
	 NEWTON = 15; ERROR = 2.4271E-03; RES = 1.3468E-05
	 NEWTON = 16; ERROR = 1.6761E-03; RES = 9.3140E-06
	 NEWTON = 17; ERROR = 1.1582E-03; RES = 6.4440E-06
	 NEWTON = 18; ERROR = 8.0040E-04; RES = 4.4561E-06
LOAD STEP = 28; TIME = 1537536.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.0481E-03
	 NEWTON = 2; ERROR = 5.7907E-01; RES = 1.4483E-03
	 NEWTON = 3; ERROR = 3.0359E-01; RES = 1.2831E-03
	 NEWTON = 4; ERROR = 1.8238E-01; RES = 8.0851E-04
	 NEWTON = 5; ERROR = 1.1821E-01; RES = 5.8300E-04
	 NEWTON = 6; ERROR = 7.8675E-02; RES = 4.0057E-04
	 NEWTON = 7; ERROR = 5.3650E-02; RES = 2.8629E-04
	 NEWTON = 8; ERROR = 3.6605E-02; RES = 2.0082E-04
	 NEWTON = 9; ERROR = 2.5376E-02; RES = 1.4235E-04
	 NEWTON = 10; ERROR = 1.7709E-02; RES = 1.0354E-04
	 NEWTON = 11; ERROR = 1.2417E-02; RES = 7.2494E-05
	 NEWTON = 12; ERROR = 8.7214E-03; RES = 5.1575E-05
	 NEWTON = 13; ERROR = 6.1413E-03; RES = 3.6435E-05
	 NEWTON = 14; ERROR = 4.3276E-03; RES = 2.5800E-05
	 NEWTON = 15; ERROR = 3.0537E-03; RES = 1.8247E-05
	 NEWTON = 16; ERROR = 2.1554E-03; RES = 1.2907E-05
	 NEWTON = 17; ERROR = 1.5224E-03; RES = 9.1284E-06
	 NEWTON = 18; ERROR = 1.0755E-03; RES = 6.4548E-06
	 NEWTON = 19; ERROR = 7.6001E-04; RES = 4.5646E-06
LOAD STEP = 29; TIME = 1592448.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.3192E-03
	 NEWTON = 2; ERROR = 5.7145E-01; RES = 1.6173E-03
	 NEWTON = 3; ERROR = 3.1238E-01; RES = 1.4906E-03
	 NEWTON = 4; ERROR = 1.8984E-01; RES = 8.8385E-04
	 NEWTON = 5; ERROR = 1.2478E-01; RES = 6.6756E-04
	 NEWTON = 6; ERROR = 8.2836E-02; RES = 4.5919E-04
	 NEWTON = 7; ERROR = 5.6979E-02; RES = 3.3324E-04
	 NEWTON = 8; ERROR = 3.9101E-02; RES = 2.3405E-04
	 NEWTON = 9; ERROR = 2.7339E-02; RES = 1.6780E-04
	 NEWTON = 10; ERROR = 1.9068E-02; RES = 1.1841E-04
	 NEWTON = 11; ERROR = 1.3432E-02; RES = 8.4479E-05
	 NEWTON = 12; ERROR = 9.4430E-03; RES = 5.9739E-05
	 NEWTON = 13; ERROR = 6.6746E-03; RES = 4.2505E-05
	 NEWTON = 14; ERROR = 4.7111E-03; RES = 3.0088E-05
	 NEWTON = 15; ERROR = 3.3353E-03; RES = 2.1375E-05
	 NEWTON = 16; ERROR = 2.3590E-03; RES = 1.5140E-05
	 NEWTON = 17; ERROR = 1.6714E-03; RES = 1.0746E-05
	 NEWTON = 18; ERROR = 1.1834E-03; RES = 7.6141E-06
	 NEWTON = 19; ERROR = 8.3874E-04; RES = 5.4016E-06
LOAD STEP = 30; TIME = 1647360.00:
	 NEWTON = 1; ERROR = 1.0000E+00; RES = 3.7058E-03
	 NEWTON = 2; ERROR = 6.2906E-01; RES = 1.7245E-03
	 NEWTON = 3; ERROR = 3.2952E-01; RES = 1.7021E-03
	 NEWTON = 4; ERROR = 2.0047E-01; RES = 9.9402E-04
	 NEWTON = 5; ERROR = 1.3036E-01; RES = 7.4563E-04
	 NEWTON = 6; ERROR = 8.8999E-02; RES = 5.1053E-04
	 NEWTON = 7; ERROR = 6.1408E-02; RES = 4.3320E-04
	 NEWTON = 8; ERROR = 4.3126E-02; RES = 2.6721E-04
	 NEWTON = 9; ERROR = 3.1083E-02; RES = 2.2538E-04
	 NEWTON = 10; ERROR = 2.2587E-02; RES = 1.5549E-04
	 NEWTON = 11; ERROR = 1.6223E-02; RES = 1.1925E-04
	 NEWTON = 12; ERROR = 1.1936E-02; RES = 8.5683E-05
	 NEWTON = 13; ERROR = 8.6429E-03; RES = 6.4495E-05
	 NEWTON = 14; ERROR = 6.3722E-03; RES = 4.6843E-05
	 NEWTON = 15; ERROR = 4.6447E-03; RES = 3.5038E-05
	 NEWTON = 16; ERROR = 3.4255E-03; RES = 2.5446E-05
	 NEWTON = 17; ERROR = 2.5060E-03; RES = 1.9028E-05
	 NEWTON = 18; ERROR = 1.8477E-03; RES = 1.3779E-05
	 NEWTON = 19; ERROR = 1.3547E-03; RES = 1.0340E-05
	 NEWTON = 20; ERROR = 9.9847E-04; RES = 7.4423E-06
```

## To do list
- Divide the contact surface into many small segments to see the effectiveness.
- Consider objective stress rate.
- Re-meshing strategy for large deformation simulation.