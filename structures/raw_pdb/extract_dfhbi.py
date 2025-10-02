from Bio.PDB import PDBParser, PDBIO, Select

class SmallMoleculeSelect(Select):
    def accept_residue(self, residue):
        return residue.get_resname() == '2ZY'

parser = PDBParser()
structure = parser.get_structure('RNA', '8K7W.pdb')
io = PDBIO()
io.set_structure(structure)
io.save('dfhbi_extracted.pdb', SmallMoleculeSelect())
