# MembraneDepthMapping
Map the depth on a PDB file. The depth, here, correspond to the distance between the upper phosphate (COM*) plane and the protein (COM).
* COM = center of mass


# Usage
python MapDepth.py -pdb my.pdb -Up [optional] -Low [optional]

- pdb: the PDB file of your peripheral membrane protein in complex with a bilayer
- Up: If your protein bound to the upperleaflet
- Low: If your protein bound to the lowerleaflet
 
output: a pdb file with the depth instaead of the B-factor. You can use Pymol or Chimera to visualize it:

![](example.png "output" )
