from dataclasses import dataclass, field


@dataclass
class PSMolecule:
	smiles: str
	HEG: int
	mp: float
	TE: float


entryOne = PSMolecule(smiles="N1N=NN=C1C1=CC=CC=C1 |c:1,3,8,10,t:6|", HEG=1, mp=136.7, TE=242.8)





print(entryOne)