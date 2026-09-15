"""Stage 1 descriptors D1-D18, exactly as frozen in PREREGISTRATION_acylation_stage1.md
(sha256 708e522fb070e2a270b8e71c23f56cca0041e6cf905cf89385be7d9a092bb28f).
Functions of sub_2_smiles only: no yield, no condition, no class label."""
import pandas as pd, numpy as np, os
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")
COOH = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

def fused_system(mol, start):
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()
             if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    sysm, frontier, seen = set(), [r for r in rings if start in r], []
    while frontier:
        r = frontier.pop()
        if r in seen: continue
        seen.append(r); sysm |= r
        for o in rings:
            if o not in seen and o & r: frontier.append(o)
    return sysm, seen

def describe(smi):
    mol = Chem.MolFromSmiles(smi)
    match = mol.GetSubstructMatches(COOH)[0]
    cidx, oc, oh = match[0], match[1], match[2]
    catom = mol.GetAtomWithIdx(cidx)
    nbrs = [n for n in catom.GetNeighbors() if n.GetIdx() not in (oc, oh)]
    a = nbrs[0]
    arom = a.GetIsAromatic()
    d = {"acid_smiles": smi}
    if arom:
        sysm, rings = fused_system(mol, a.GetIdx())
        hets = [mol.GetAtomWithIdx(i).GetSymbol() for i in sysm if mol.GetAtomWithIdx(i).GetSymbol() != "C"]
        own = [r for r in rings if a.GetIdx() in r][0]
        own_hets = [mol.GetAtomWithIdx(i).GetSymbol() for i in own if mol.GetAtomWithIdx(i).GetSymbol() != "C"]
        frag = Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(sysm), canonical=True)
        d["D1_ring_system_state"] = "A" if hets else "B"
        d["D4_attached_ring_hetero"] = int(bool(own_hets))
        d["D6_hetero_N"] = int("N" in hets); d["D6_hetero_O"] = int("O" in hets); d["D6_hetero_S"] = int("S" in hets)
        d["D7_n_ring_heteroatoms"] = len(hets)
        d["D8_ring_system_frag"] = frag
        d["D14_n_ortho_substituents"] = sum(
            1 for n in a.GetNeighbors() if n.GetIdx() in sysm
            for m in n.GetNeighbors() if m.GetIdx() not in sysm and m.GetIdx() != cidx)
        d["D10_hetero_adjacent"] = int(any(
            n.GetSymbol() != "C" and n.GetIsAromatic() for n in a.GetNeighbors()))
    else:
        d["D1_ring_system_state"] = "C"
        d["D4_attached_ring_hetero"] = 0
        d["D6_hetero_N"] = d["D6_hetero_O"] = d["D6_hetero_S"] = 0
        d["D7_n_ring_heteroatoms"] = 0
        d["D8_ring_system_frag"] = "aliphatic"
        d["D14_n_ortho_substituents"] = sum(
            1 for n in a.GetNeighbors() if n.GetIdx() != cidx
            for m in n.GetNeighbors() if m.GetIdx() not in (a.GetIdx(), cidx))
        d["D10_hetero_adjacent"] = 0
    d["D2_aromatic_attachment"] = int(arom)
    d["D3_hetero_vs_rest"] = int(d["D1_ring_system_state"] == "A")
    d["D5_molecule_has_arom_hetero"] = int(any(
        at.GetIsAromatic() and at.GetSymbol() != "C" for at in mol.GetAtoms()))
    # D9 shortest bond path from carboxyl C to nearest AROMATIC ring heteroatom
    dm = Chem.GetDistanceMatrix(mol)
    tgt = [at.GetIdx() for at in mol.GetAtoms() if at.GetIsAromatic() and at.GetSymbol() != "C" and at.IsInRing()]
    d["D9_dist_to_ring_hetero"] = int(min(dm[cidx][t] for t in tgt)) if tgt else -1
    mh = Chem.AddHs(mol)
    from rdkit.Chem import AllChem
    AllChem.ComputeGasteigerCharges(mh)
    d["D11_q_hydroxyl_O"] = float(mh.GetAtomWithIdx(oh).GetDoubleProp("_GasteigerCharge"))
    d["D12_q_carbonyl_C"] = float(mh.GetAtomWithIdx(cidx).GetDoubleProp("_GasteigerCharge"))
    d["D13_alpha_heavy_degree"] = a.GetDegree()
    d["D15_n_rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    d["D16_MolWt"] = Descriptors.MolWt(mol)
    d["D17_TPSA"] = rdMolDescriptors.CalcTPSA(mol)
    d["D18_MolLogP"] = Crippen.MolLogP(mol)
    sc = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    d["murcko"] = sc if sc else "ACYCLIC"
    return d

# --- ring-system naming, frozen with the descriptor list ---
_NAME = {"aliphatic": "aliphatic", "c1ccccc1": "benzene", "c1ccc2ccccc2c1": "naphthalene",
         "c1ccncc1": "pyridine", "c1cnccn1": "pyrazine", "c1cncnc1": "pyrimidine",
         "c1ccsc1": "thiophene", "c1ccoc1": "furan", "c1ccc2ccccc2n1": "quinoline",
         "c1ccc2nccnc2c1": "quinoxaline", "c1ccc2[nH]ccc2c1": "indole",
         "c1ccc2occc2c1": "benzofuran", "c1ccc2sccc2c1": "benzothiophene",
         "c1cn[nH]c1": "pyrazole", "c1c[nH]cn1": "imidazole", "c1cnc[nH]1": "imidazole"}


def ring_system_name(frag):
    if frag == "aliphatic":
        return "aliphatic"
    try:
        c = Chem.CanonSmiles(frag)
    except Exception:
        c = frag
    return _NAME.get(c, _NAME.get(frag, c))


def descriptor_table(smiles_list):
    """The frozen D1-D18 block plus the Murcko scaffold, one row per acid."""
    import pandas as pd
    df = pd.DataFrame([describe(s) for s in smiles_list])
    df["D8_ring_system_class"] = df.D8_ring_system_frag.map(ring_system_name)
    return df
