import sys
from Bio.PDB import PDBParser, NeighborSearch

def get_interface_residues(pdb_file, target_chain_id='A', cutoff=4.0):
    """
    Identifies residues in target_chain_id that are within 'cutoff' Angstroms
    of any atom in any other chain.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)
    model = structure[0]

    target_atoms = []
    other_atoms = []

    for chain in model:
        if chain.id == target_chain_id:
            for residue in chain:
                if residue.id[0] == ' ': 
                    target_atoms.extend(residue.get_atoms())
        else:
            for residue in chain:
                if residue.id[0] == ' ':
                    other_atoms.extend(residue.get_atoms())

    if not other_atoms:
        print("Error: No other chains found. Is this already a monomer?")
        return []

    ns = NeighborSearch(other_atoms)
    interface_residues = set()

    print(f"Scanning Chain {target_chain_id} against all other chains...")
    print(f"Distance Cutoff: {cutoff} Angstroms")

    for atom in target_atoms:
        neighbors = ns.search(atom.coord, cutoff)
        if neighbors:
            res_id = atom.get_parent().id[1]
            interface_residues.add(res_id)

    sorted_interface = sorted(list(interface_residues))
    
    print(f"\nFound {len(sorted_interface)} interface residues.")
    print("-" * 40)
    print("COPY THIS LIST INTO YOUR SCANNING SCRIPT:")
    print(f"INTERFACE_RESIDUES = {sorted_interface}")
    print("-" * 40)

    return sorted_interface

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python find_interface_residues.py <tetramer_pdb>")
        sys.exit(1)
    
    get_interface_residues(sys.argv[1])
