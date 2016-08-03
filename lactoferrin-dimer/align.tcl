# Molecule 0 is the reference structure (CG)
# Molecule 1 is the structure to overlay. (all atom)
set sel2 [atomselect 0 "resname GLU and index >= 697"]
set sel1 [atomselect 2 "resname GLU and name CB"]
set sel3 [atomselect 2 "all"]

set transformation_matrix [measure fit $sel1 $sel2]

$sel3 move $transformation_matrix
