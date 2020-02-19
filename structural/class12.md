Visualizing Protein Structure
================
Vikas Gubbala
2/13/2020

## Structural Analysis of HIV-1 Protease

HIV-1 protease is essential for the formation of HIV proteins.
Inhibition of this protein may aid in the suppress viral replication. In
this analysis, we will look at the PDB file 1HSG, which contains the
protease bound to the drug indinavir.

PDB files contain lists of atoms along with coordinates. Sometimes,
additional information like secondary structure is added.

``` r
library(bio3d)
library(ggplot2)
library(ggrepel)
```

## Working with Single PDB files.

The 1HSG file is read from the online database.

``` r
#accesses online pdb file
pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
# summary of pdb object
pdb 
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

We can select for certain attributes as well.

``` r
#attributes of PDB object
attributes(pdb) 
```

    ## $names
    ## [1] "atom"   "xyz"    "seqres" "helix"  "sheet"  "calpha" "remark" "call"  
    ## 
    ## $class
    ## [1] "pdb" "sse"

``` r
#selects an attribute of the object
head(pdb$atom) 
```

    ##   type eleno elety  alt resid chain resno insert      x      y     z o     b
    ## 1 ATOM     1     N <NA>   PRO     A     1   <NA> 29.361 39.686 5.862 1 38.10
    ## 2 ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1 40.62
    ## 3 ATOM     3     C <NA>   PRO     A     1   <NA> 29.760 38.071 4.022 1 42.64
    ## 4 ATOM     4     O <NA>   PRO     A     1   <NA> 28.600 38.302 3.676 1 43.40
    ## 5 ATOM     5    CB <NA>   PRO     A     1   <NA> 30.508 37.541 6.342 1 37.87
    ## 6 ATOM     6    CG <NA>   PRO     A     1   <NA> 29.296 37.591 7.162 1 38.40
    ##   segid elesy charge
    ## 1  <NA>     N   <NA>
    ## 2  <NA>     C   <NA>
    ## 3  <NA>     C   <NA>
    ## 4  <NA>     O   <NA>
    ## 5  <NA>     C   <NA>
    ## 6  <NA>     C   <NA>

``` r
#subset atoms
pdb$atom[1:5, c("eleno", "elety")]
```

    ##   eleno elety
    ## 1     1     N
    ## 2     2    CA
    ## 3     3     C
    ## 4     4     O
    ## 5     5    CB

Protein Visualization is done in the software VMD.

``` r
## Do this!!
```

protein \<- atom.select(pdb, “protein”, value = TRUE) ligand \<-
atom.select(pdb, “ligand”, value = TRUE)
