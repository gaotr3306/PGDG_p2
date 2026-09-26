# PGDG_p2

Source codes accompanying the paper:

**A Particle-Grid Discontinuous Galerkin Method, Part 2:
Strongly Compressible Multiphase Flows**

## Examples

This repository contains five standalone numerical examples:

1. 3D Air/He shock-bubble interaction (cylindrical bubble)
2. 3D Air/He shock-bubble interaction (spherical bubble)
3. Underwater explosion
4. Triple-point problem (two-phase)
5. Triple-point problem (three-phase)

## Code Organization

The five numerical examples are currently provided as standalone programs.
They share the same PG-DG methodology but have not yet been organized into
a unified general-purpose solver.

The `3D_generalized` directory provides a relatively generalized implementation
of the 3D PG-DG framework. Users who wish to extend the present code to other
problems are recommended to start from this version.

In `3D_generalized`, the function `Split_particle_fast` includes specific
treatments for inlet and outlet particles. Users should modify this part
accordingly for problems without inlet or outlet boundaries.

## Dependencies

The following libraries were used in the present implementation:

- oneTBB 2022.0.0
- MPICH 3.4.1
- Boost 1.80.0
- Eigen 3.4.0

Other compatible versions may also work but have not been systematically tested.

## Compilation and Usage

Each numerical example is provided as a standalone program. Please refer to
the `README.md` file in each example directory for case-specific compilation,
configuration, and execution instructions.

## Citation

If you use this code in your research, please cite the corresponding paper:

**T. Gao and D. Hyde, "A Particle-Grid Discontinuous Galerkin Method,
Part 2: Strongly Compressible Multiphase Flows."**

## License

[License information to be added.]