# GEOKEYFEM_HM
The numerical simulation code of Geoinvention group (PI: Prof. Zhenyu YIN) of PolyU, mainly developed by Qi ZHANG.

## Functionality
The code simulates a contact problem between a rigid rectangular block with a Mohr-Coulomb soil by using the penalty method (**small deformation**). The deformation equation $\sigma_{ij,j} = 0$ is discretized  by using the [Smoothed Finite Element Method](https://www.taylorfrancis.com/books/mono/10.1201/EBK1439820278/smoothed-finite-element-methods-liu-nguyen-trung). The direct nodal integration on the smoothing domain is also modified to **stabilized conforming nodal integration (SCNI)**. This is reflected in this [file](./assemble_stab.m).


Although the current example only considers deformation field, the code is designed for hydromechanical coupling analysis (poromechanics). Therefore, the [pore pressure stabilization technique](https://doi.org/10.1016/j.cma.2008.05.015) is included. In addition, for simplicity, the biot coefficient $b$ is set to be `1` and Biot modulus $M$ is `infinite`. It should be not too difficult to modify our code for a more general case by using relevant pre-defined matrices such as the `Mass matrix` $\int_{\Omega} N^T N {\rm d}\ V$.


This is the [main file](./main_rigid_contact_prob.m). You can run it directly. In the main file, the most time-consuming part is the assembly process, whose function is given in this [file](./assemble_system.m). The [constitutive file](./MohrCoulomb_UMAT.m) is also called in the assembly process, as shown in the following code block:
```
[stress_new(:, ino), SDV_new(:, ino), cto] = MohrCoulomb_UMAT(0, Props(ino,:), stress(:, ino), strain_ino_new - strain_ino_old, SDV(:, ino));
```

The [assign_tractionBC2](./assign_tractionBC2.m) is not used in this contact problem, while it is designed to calculate the equivalent nodal force vector in FEM. The user only need to define the `traction_f` as a function of both location *x* and time *t*.


## Output
