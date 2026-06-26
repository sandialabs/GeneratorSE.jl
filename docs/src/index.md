# GeneratorSE.jl

`GeneratorSE` provides analytical sizing and performance routines for permanent
magnet synchronous generators. The package is a Julia translation of the NREL
GeneratorSE-style models, with exported functions for radial-flux, axial-flux,
Halbach axial-flux, dynamic PMSG behavior, and generator costing.

In the SIRENOpt ontology, GeneratorSE is used as the electromechanical boundary
between rotor shaft power and electrical power. The stable integration quantities
are torque, speed, generated power, efficiency/losses, mass, and cost.
