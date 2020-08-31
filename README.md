# cannabinoid1_scra

Steps to reproduce:
- Run `align_drugs.py`. This uses rdkit to generate 1000 conformers of each SCRA ligand, then aligns them to MDMB-Fubinaca from the 6N4B crystal structure. That way, the starting coordinates closely match the crystallized ligand.
- Run `parameterize_drugs.py`. This loads the PDB files for each ligand and parameterizes it with OpenForceField.
- Run `do_equilibrations.sh`. This runs python script `equilibrate.py` for the ligands of interest. The equilibration consists of applying restraints to the protein coordinates and gradually relaxing. Relaxation of restraints means every 12ps `restraint = 0.9*restraint`
- Run `do_productions.sh`. This runs python script `production.py` which simulates for 90ns using OpenMM.
