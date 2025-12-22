import math

pdb_file = "1GV3_clean_Repair_11.pdb"

def get_distance(coord1, coord2):
    return math.sqrt(
        (coord1[0] - coord2[0])**2 +
        (coord1[1] - coord2[1])**2 +
        (coord1[2] - coord2[2])**2
    )

print(f"--- Reading {pdb_file} ---")

atoms = []
try:
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                chain = line[21]
                res_seq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                if atom_name.startswith("N") or atom_name.startswith("O"):
                    atoms.append({
                        "type": "N" if atom_name.startswith("N") else "O",
                        "chain": chain,
                        "res": res_seq,
                        "coords": (x, y, z)
                    })

    print(f"Found {len(atoms)} N/O atoms. Calculating interactions...")

    hbond_count = 0
    
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            a1 = atoms[i]
            a2 = atoms[j]
            
            if a1['type'] == a2['type']: continue
            
            if a1['chain'] == a2['chain'] and a1['res'] == a2['res']: continue
            
            dist = get_distance(a1['coords'], a2['coords'])
            
            if 2.4 < dist < 3.5:
                hbond_count += 1

    print("-" * 30)
    print(f"Estimated H-Bonds (Heavy Atom Proxy): {hbond_count}")
    print("-" * 30)

except FileNotFoundError:
    print(f"Error: Could not find {pdb_file}")
