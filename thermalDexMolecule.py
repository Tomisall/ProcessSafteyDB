from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles
from dataclasses import dataclass, field
from io import BytesIO
import re

highEnergyGroups = ".\\_core\\HighEnergyGroups.csv"
expEnergyGroups = ".\\_core\\ExplosiveGroups.csv"

@dataclass
class thermalDexMolecule:
    # Class for keeping track of Molecules in ThermalDex

    # User Provided Values
    SMILES: str
    name: str = ''
    mp: float = None
    mpEnd: float = None
    Q_dsc: float = None
    onsetT: float = None
    initT: float = None
    proj: str = ''

    # Calculated Values
    mol: any = None
    MW: float = None
    molForm: str = None
    eleComp: str = None
    HEG: int = None
    HEG_list: list = field(default_factory=list)
    EFG: int = None
    EFG_list: list = field(default_factory=list)
    RoS_val: float = None
    RoS_des: str = None
    OB_val: float = None
    OB_des: str = None
    IS_val: float = None
    IS_des: str = None
    EP_val: float = None
    EP_des: str = None
    Td24: float = None
    oreoSmallScale_val: int = None
    oreoSmallScale_des: str = None
    oreoTensScale_val: int = None
    oreoTensScale_des: str = None
    oreoHundsScale_val: int = None
    oreoHundsScale_des: str = None
    oreoLargeScale_val: int = None
    oreoLargeScale_des: str = None

    def genMol(self):
        RDmol = MolFromSmiles(self.SMILES)
        self.mol = RDmol
        return RDmol

    def molToQPixmap(self):
        # Generate a molecular drawing as a PNG image
        img = Draw.MolToImage(self.mol)

        # Convert the image to a byte array
        byte_array = BytesIO()
        img.save(byte_array, format='PNG')

        # Convert the byte array to a QPixmap and display it
        pixmap = QPixmap()
        pixmap.loadFromData(byte_array.getvalue())
        return pixmap

    def mwFromMol(self):
        cmpdMW = Descriptors.MolWt(self.mol)
        self.MW = cmpdMW

    def HEGFromMol(self):
        fullMatch = 0
        with open(highEnergyGroups, "r") as HEGroups: 
            for line in HEGroups:
                HeSubstructure = MolFromSmiles(line)
                fullmatchList = Mol.GetSubstructMatches(self.mol, HeSubstructure)
                if len(fullmatchList) > 0:
                    print('High Energy Group Found: ' + line[:-1])
                    self.HEG_list += [line[:-1]]
                fullMatch += len(fullmatchList)

        self.HEG = fullMatch

    def EFGFromMol(self):
        expMatch= 0
        with open(expEnergyGroups, "r") as expGroups: 
            for line in expGroups:
                expSubstructure = MolFromSmiles(line)
                expmatchList = Mol.GetSubstructMatches(self.mol, expSubstructure)
                if len(expmatchList) > 0:
                    print('Explosive Group Found: ' + line[:-1])
                    self.EFG_list += [line[:-1]]
                expMatch += len(expmatchList)

        self.EFG = expMatch   

    def eleCompFromMol(self):
        self.molForm = rdMolDescriptors.CalcMolFormula(self.mol)
        eleComp = ""
        eleCompList = []
        eleList = []
        niceList = []
        atomPattern = r"Cl\d*|H\d*|O\d*|N\d*|Si\d*|S\d*|F\d*|Cs\d*|Br\d*|I\d*|B\d*|Al\d*|Na\d*|K\d*|Mg\d*|Zn\d*|Ti\d*|Pd\d*|C\d*"
        match = re.findall(atomPattern, self.molForm, re.I)
        if match:
            items = match
        for ele in match:
            ment = re.findall(r"([a-z]+)([0-9]+)?", ele, re.I)
            matchedEleComp = ment[0]
            eleList += [matchedEleComp]
            
        for compostion in eleList:
            eleCompList += compostion[::-1]
            niceList += [('').join(compostion[::-1])]
   
        self.eleComp = (', ').join(niceList)
           


nicatinamide = thermalDexMolecule("ClC(=O)C1=CC=NC=C1 |c:5,7,t:3|")
nicatinamide.genMol()
nicatinamide.mwFromMol()
nicatinamide.HEGFromMol()
nicatinamide.EFGFromMol()
nicatinamide.eleCompFromMol()
print('\n')
print(nicatinamide)