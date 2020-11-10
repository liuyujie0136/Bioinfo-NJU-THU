### Using bio3d to do biological structure analysis of protein structure and sequence data

library(bio3d)

## Try some demos:
demo(package="bio3d") # all are listed below
demo(pdb) # PDB Reading, Manipulation, Searching and Alignment
demo(pca) # Principal Component Analysis
demo(md)  # Molecular Dynamics Trajectory Analysis
demo(nma) # Normal Mode Analysis


## Read a PDB file
pdb <- read.pdb("4q21")  # From local file, use: "xxx.pdb"


## Print a brief composition summary
pdb
attributes(pdb)


## Examine the storage format (or internal *str*ucture)
str(pdb)


## Print data for the first four atom
pdb$atom[1:4,]


## PLOT: Residue temperature factors for PDB file with secondary structure element (SSE) annotation in marginal regions
plot.bio3d(pdb$atom$b[pdb$calpha], sse=pdb, typ="l", ylab="B-factor")


## Print some coordinate data
head(pdb$atom[, c("x","y","z")])

# OR coordinates as a numeric vector
pdb$xyz
pdb$xyz[,1:3]


## Print C-alpha coordinates
head(pdb$atom[pdb$calpha, c("resid","elety","x","y","z")])

# OR use 'atom.select' function
ca.inds <- atom.select(pdb, elety="CA") # returns 'indices' (row numbers)
head(pdb$atom[ca.inds$atom,])
head(pdb$xyz[,ca.inds$xyz])


## More examples about atom.select()
# string (calpha, cbeta, backbone, sidechain, protein, nucleic, ligand, water, h, or noh); chain; resno; inverse...
nuc.inds <- atom.select(pdb,"nucleic")
a.inds <- atom.select(pdb, chain="A")
cab.inds <- atom.select(pdb, elety=c("CA","CB"), chain="A", resno=10:20)
nowat.inds <- atom.select(pdb, "water", inverse=TRUE)

# Combining selections
at.conb.sele <- atom.select(pdb, "protein", resid="GDP", operator="OR")
a.inds <- atom.select(pdb, string="protein")
b.inds <- atom.select(pdb, resid="GDP")
conb.sele <- combine.select(a.inds, b.inds, operator="OR")


## Extract amino acid sequence in 3-letter and 1-letter forms
# use atom.select()
ca.inds <- atom.select(pdb, "calpha")
aa3 <- pdb$atom$resid[ca.inds$atom]
head(aa3)
aa1 <- aa321(aa3)  # 3-letters to 1-letter
head(aa1)

# OR use pdbseq()
aa1 <- pdbseq(pdb)
aa3 <- aa123(aa1)


## Output a backbone only PDB file
b.inds <- atom.select(pdb, "back")
backpdb <- trim.pdb(pdb, b.inds)
write.pdb(backpdb, file="backbone_only.pdb")

# Selection statements can be passed directly to trim.pdb()
backpdb <- trim.pdb(pdb, "backbone")
# The 'value=TRUE' option of atom.select() will return a PDB object
backpdb <- atom.select(pdb, "backbone", value=TRUE)


## Example of Manipulate a PDB object
# select chains A, E and F
inds <- atom.select(pdb, chain=c("A", "E", "F"))

# trim PDB to selection
pdb2 <- trim.pdb(pdb, inds)

# assign new chain identifiers
pdb2$atom$chain[ pdb2$atom$chain=="E" ] <- "B"
pdb2$atom$chain[ pdb2$atom$chain=="F" ] <- "C"

# re-number chain B and C
pdb2$atom$resno[ pdb2$atom$chain=="B" ] <- pdb2$atom$resno[ pdb2$atom$chain=="B" ] - 156
pdb2$atom$resno[ pdb2$atom$chain=="C" ] <- pdb2$atom$resno[ pdb2$atom$chain=="C" ] - 156

# assign the GDP residue a residue number of 500
pdb2$atom$resno[ pdb2$atom$resid=="GDP" ] <- 500

# use chain D for the GDP residue
pdb2$atom$chain[ pdb2$atom$resid=="GDP" ] <- "D"

# Center, to the coordinate origin, and orient, by principal axes,
# the coordinates of a given PDB structure or xyz vector.
xyz <- orient.pdb(pdb2)

# write the new pdb object to file
write.pdb(pdb2, xyz=xyz, file="4LHY_AEF-oriented.pdb")


## Coordinate superposition and structural alignment (need MUSCLE!!)
# Align and superpose two or more structures
pdbs <- pdbaln(c("4q21", "521p"), fit=TRUE)
pdbs


## Binding site identification
# read G-protein structure
pdb <- read.pdb("4q21")
bs <- binding.site(pdb)

# residue names of identified binding site between protein and ligand
bs$resnames


## Example of determination of residues at the binding interface between two proteins
b <- read.pdb("4lhy")

# atom selection
a.inds <- atom.select(b, chain="A")
b.inds <- atom.select(b, chain=c("E", "F"))

# identify interface residues
bs <- binding.site(b, a.inds=a.inds, b.inds=b.inds)

# use b-factor column to store interface in PDB file
b$atom$b[ bs$inds$atom ] <- 1
b$atom$b[ -bs$inds$atom ] <- 0

# write to file
write.pdb(b, file="4LHY-interface.pdb")


## Reading multi-model PDB files
# Read multi-model PDB file
pdb.multi <- read.pdb("1d1d", multi=TRUE)

# The xyz component contains 20 frames
pdb.multi$xyz

# Select a subset of the protein
ca.inds <- atom.select(pdb.multi, "calpha")

# Access C-alpha coordinates of the first 5 models
pdb.multi$xyz[1:5, ca.inds$xyz]


## Identification of dynamic domains
# Domain analysis
gs  <- geostas(pdb.multi)

# Fit all frames to the 'first' domain
domain.inds <- gs$inds[[1]]
xyz <- pdbfit(pdb.multi, inds=domain.inds)

# write fitted coordinates
write.pdb(pdb.multi, xyz=xyz, chain=gs$atomgrps, file="1d1d_fit-domain1.pdb")

# plot geostas results
plot(gs, contour=FALSE)


## Invariant core identification
# Invariant core
core <- core.find(pdb.multi)

# fit to core region
xyz <- pdbfit(pdb.multi, inds=core)

# write fitted coordinates
write.pdb(pdb.multi, xyz=xyz, file="1d1d_fit-core.pdb")


## Constructing biological units
pdb$remark$biomat
biounit(pdb)


## Working with multiple PDB files - basic
# Download some example PDB files
ids <- c("1TND_B","1AGR_A","1FQJ_A","1TAG_A","1GG2_A","1KJY_A")
raw.files <- get.pdb(ids)

# Extract and align the chains we are interested in
files <- pdbsplit(raw.files, ids)
pdbs <- pdbaln(files)

# Calculate sequence identity
pdbs$id <- basename.pdb(pdbs$id)
seqidentity(pdbs)

# Calculate RMSD
rmsd(pdbs, fit=TRUE)

# Quick PCA
pc <- pca(pdbfit(pdbs), rm.gaps=TRUE)
plot(pc)

# Quick NMA of all structures
modes <- nma(pdbs)
plot(modes, pdbs, spread=TRUE)

# Extracting SSE (secondary structure element)
pdbs$sse

