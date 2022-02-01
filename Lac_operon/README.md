The Lactose operon (denoted gop or Gop) found in E. coli consists: of operon-associated promoter and operator sites; and of three genes (denoted lacAYZ) encoding the three proteins (denoted LacAYZ) responsible for the import and first steps of lactose metabolism. Gop expression can be represented by the equation: Gop -> Gop + LacAYZ 

Gop is downregulated by the binding of repressor protein LacI, which is continuously expressed by the lacI gene. LacI also binds to allolactose (denoted alac) which is a product of lactose metabolism in the cell. In presence of lactose, alac is synthesised and LacI binds alac in preference to gop, such that the cell is able to express LacAYZ from gop and metabolise the lactose. 

MATLAB code to model the switching behaviour of the lactose operon.
Set rate constants (k1-k4) in the lac_rep_model function and the initial species concentrations (x1-x5) in the ode_45_example_simulation files.
Run the simulation from ode_45_example_simulation to produce a representative graph of species dynamics at set conditions.
