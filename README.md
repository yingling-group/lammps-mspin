# LAMMPS-MSPIN
A LAMMPS plugin `MSPIN` for atomistic molecular dynamics simulations of magnetic nanoparticles.

Developed for simulation method published in:
>A.U. Mahmood and Y.G. Yingling. *All-Atom Simulation Method for Zeeman Alignment
and Dipolar Assembly of Magnetic Nanoparticles.* (2021)

Please cite paper above if used.

Author:   
Akhlak Mahmood   
Yingling Group  
North Carolina State University

## Branches
- If you are interested in installing the `MSPIN` plugin to your local machine,
please see the instructions on the [main branch](https://github.com/yingling-group/lammps-mspin "main branch").

- If you are interested in the source code of the plugin, please see the 
[plugin branch](https://github.com/yingling-group/lammps-mspin/tree/plugin "plugin branch").


# INSTALLATION

1. Clone the main branch to your local machine.
    ```
    git clone --branch main https://github.com/yingling-group/lammps-mspin.git
    ```

2. CD to the cloned directory and create a build directory.
    ```
    cd lammps-mspin
    mkdir build
    ```
3. Run CMake with the `user.cmake` preset from the build directory.
    ```
    cd build
    cmake3 -C ../cmake/presets/user.cmake ../cmake
    ```

4. Build.
    ```
    make -j4
    ```

# USAGE

General steps to follow:
1. Identify two atoms in the MNP core to set as 'MS' atoms.
2. Calculate the distance between the two atoms, required "magnetic charge", Qm.
3. Generate a LAMMPS data file using the newly defined `qmag` atom style with the calculated Qm values.
Make sure that one of the "MS" atoms has a positive Qm value, and the other has a negative Qm value.
2. Create a LAMMPS input file and define magnetic interactions using the commands described below.
3. Run the simulation.

## Data File for MNP
A new atom style `qmag` has been introduced in the plugin which is the LAMMPS `full` atom style
with an additional column/field for "magnetic charge" Qm values.
An example of a data file can be found in the `examples/mspin` directory.

### "Magnetic Charge" Calculation
As an example, we will consider the magnetic moment of a Fe3O4 unit cell, which is `32 muB`,
where `muB` is the Bohr Magneton with a value,

				muB 	= 9.274E-24 Joules/Tesla
    					= 9.274E-24 Ampere Meter^2
	    				= 5.7884 e / fs Angstrom^2

Also, the lattice constant of Fe3O4 is `8.39 Angstrom` and
magnetization is defined as the dipole moment per unit volume.
So, we can find the volume of a Fe3O4 unit cell `= (8.39 Angstrom)^3 = 5.905E-28 Meter^3`.
    
Then the magnetization,

				M	= 32 * muB /  5.905E-28 Meter^3
					= 502.57E+3 Ampere/Meter

To use this magnetization value, we need to convert it into LAMMPS real units.

			1 Ampere 	= 1 Coulomb / Second
    					= 6241.5 e / fs

Where, `1fs = 1E-15 seconds`, `e` is the electronic charge, `1 Coulomb = 6.2415E+18 e`.

> Note: e/fs is Ampere in LAMMPS real units.

Using this, in real units,

				M	= 0.3137 e/fs / Angstrom

Therefore, a 7nm nanoparticle will have dipole moment

				m	= M * Volume of the particle
					= 0.3137 * 4/3 * PI * 35^3 e/fs Angstrom^2
					= 56.339E+3 e/fs Angstrom^2
					= 56.339 Ke/fs Angstrom^2

If we choose two MS atoms at far edges of a nanoparticle,
then we can assume the distance between them is same as the diameter of the nanoparticle,

				d 	= 70 Angstrom

Then for each of these atoms,

				Qm 	= 56.339 / 70 Ke/fs Angstrom
					= 0.8048 Ke/fs Angstrom

Then we need to set +0.8048 and -0.8048 as +/-Qm respectively for the MS atoms of the configuration.


## Available LAMMPS Commands
To define a group for the rigid magnetic molecules. The following command defines the 
group `feo` for eight MNPs with IDs 1 to 8.

    group feo molecule 1:8

Set the MSPIN fix with the desired temperature.

    fix <fix_id> feo rigid/mspin molecule temp 300 300 10 

By default no dipole-dipole interaction is implemented.
You need to set a cutoff for the interaction using the `dpcut` keyword to turn it on.
For example, for a 64 angstroms cutoff,

    fix <fix_id> feo rigid/mspin molecule temp 300 300 10 dpcut 64

To apply an external magnetic field, use the `bfield` keyword.

    fix <fix_id> feo rigid/mspin molecule temp 300 300 10 bfield 5.0 0.0 0.0

This sets a *non-uniform* magnetic field of 5 Tesla in the x direction.
Note that, the default external B field is non-uniform.
To set a uniform magnetic field, use the `uniform` keyword.

    fix <fix_id> feo rigid/mspin molecule temp 300 300 10 bfield 5.0 0.0 0.0 uniform dpcut 64


## Thermo Outputs
Turn on the energy contributions using the `fix_modify` command.

    fix_modify <fix_id> energy yes

The Zeeman energy due to external magnetic field and the dipole-dipole interaction energy
can be separately printed. First define a new mspin type compute.

    compute magE feo mspin <fix_id>

Then use the compute variables to output them using a `thermo_style`.

    thermo_style step pe ke etot temp c_magE[1] c_magE[2]

Here, the first item `c_magE[1]` will output the total dipole-dipole interaction energy
and the second item, `c_magE[2]` will print the total Zeeman energy due to the external field.

If you have multiple mspin molecules, you can compute their center to center distance using
the `rigid/mspin/distance` compute during simulation. Example,

    compute dist feo rigid/mspin/distance 1 2

will calculate the distace between the molecules 1 and 2.
Output the value using `c_dist[1]` keyword using a `thermo_style`.

Since `MSPIN` extends the RIGID package, you can also output the translational and rotational energies
of the rigid bodies.

    compute 	    te	feo ke/rigid <fix_id>
    compute 	    re	feo erotate/rigid <fix_id>
    thermo_style 	custom step ecoul evdwl pe ke etotal temp c_te c_re

> Note: These `compute` and `thermo_style` should be defined only once
for all subsequent runs. LAMMPS doesn't permit defining compute with the
same ID more than once.

## Dipolar Interaction Scaling Factor `alpha`
Due to particles distribution and geometry, a second ‘demagnetizing factor’ comes up.
Allia et. al. Physical Review B 2001, 64 (14), 144420 mentions this term could be as much as a hundred. 
Note that this factor doesn’t affect the zeeman interaction.
When the particles are in liquid and free to move, this becomes a ‘magnetizing effect’.
See Sanchez and Raap et. al. Physical Review B 2017, 95(134421) for more info.

We can specify the value of `alpha` as a keyword while defining the MSPIN fix, (default is 1.0).

    fix <fix_id> feo rigid/mspin molecule temp 300 300 10 bfield 0.5 0.0 0.0 uniform dpcut 64 alpha 5.0

> Note: In the current implementation, `alpha` is used to scale the `mu_0/4pi` term which does not exist in the Zeeman equations.

## Dipole Moment Scaling Factor `beta`
The scaling factor `beta` is used to scale up the dipole moment directly
so it scales both Zeeman and dipolar interactions.
Since both the dipolar interaction energy and the dipole moment has
a relation with the volume of the particle, the effect of `beta` is same for both.

```
fix <fix_id> feo rigid/mspin molecule temp 300 300 10 bfield 0.5 0.0 0.0 uniform dpcut 64 alpha 5.0 beta 8.0
```

**Note**: `alpha` affects the dipolar interactions only, `beta` affects both Zeeman and dipolar interactions.

## Currently Known Limitations
- Constant pressure (NPT) simulation has not been tested.
- Only MPI parallelization is supported for the magnetic interactions.
K-space and pairwise interactions can be calculated using GPUs.
- Only real units are supported and/or tested.
- Timestep needs to be `1.0 fs`
