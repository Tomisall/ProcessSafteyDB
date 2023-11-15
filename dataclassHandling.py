from dataclasses import dataclass, field
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import csv

@dataclass
class PSMolecule:
	smiles: str = None
	HEG: int = None
	mp: float = None
	MW: float = None
	TE: float = None

	def toIterable(self):
		return iter(
			[
				self.smiles,
				self.HEG,
				self.mp,
				self.MW,
				self.TE,
			]
	)

	def toHeader(self):
		return [
			"Molecule",
			"Number of High Energy Groups",
			"Melting point",
			"Molecular Weight",
			"Thermal Event",
			]





entryOne = PSMolecule(smiles="N1N=NN=C1C1=CC=CC=C1 |c:1,3,8,10,t:6|", HEG=1, mp=136.7, TE=242.8)
RDMol = Chem.MolFromSmiles(entryOne.smiles)
cmpdMW = Descriptors.ExactMolWt(RDMol)

print(RDMol)

entryOne = PSMolecule(smiles="N1N=NN=C1C1=CC=CC=C1 |c:1,3,8,10,t:6|", HEG=1, mp=136.7, MW=cmpdMW, TE=242.8)
print(entryOne)

listforDB = [entryOne]

def writetoCSV(PSDBList: list):
	with open("PSdatabase.csv", "a") as DB:
		writer = csv.writer(DB)
		writer.writerow(PSMolecule().toHeader())
		for item in PSDBList:
			writer.writerow(item.toIterable())

writetoCSV(listforDB)

Draw.ShowMol(RDMol)
