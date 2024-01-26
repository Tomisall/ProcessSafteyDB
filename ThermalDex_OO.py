import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTableWidget, QTableWidgetItem, QTabWidget, QGraphicsPixmapItem,  QMessageBox, QComboBox #, QTableView, QToolTip
from PyQt5.QtGui import QPixmap, QColor, QIcon, QRegExpValidator #QValidator #, QCursor
from PyQt5.QtCore import Qt, QRegExp, pyqtSignal
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles
from io import BytesIO
from numpy import log10
from pubchempy import get_compounds
from dataclasses import dataclass, field
import pandas as pd
import re
import pyperclip

versionNumber = "0.5.3"

try:
    import pyi_splash
    pyi_splash.close()
except:
    pass

def importConfig():
    conf = open('.\\_core\\ThermalDex.config', 'r')
    confCounter = 0
    for line in conf:
        #print(confCounter)
        if confCounter == 4:
           defaultDB = line.strip("\n")
           confCounter += 1
        elif confCounter == 8:
           highEnergyGroups = line.strip("\n")
           confCounter += 1
        elif confCounter == 12:
           expEnergyGroups = line.strip("\n")
           confCounter += 1
        else:
           confCounter += 1

    #print(defaultDB)
    return defaultDB, highEnergyGroups, expEnergyGroups

@dataclass
class thermalDexMolecule:
    # Class for keeping track of Molecules in ThermalDex

    # User Provided Values
    SMILES: str
    name: str = ''
    mp: float = None
    mpEnd: float = None
    Q_dsc: float = None
    Qunits: str = 'J g<sup>-1</sup>'
    onsetT: float = None
    initT: float = None
    proj: str = ''

    # Calculated Values
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
    oreo: int = 0
    oreoSmallScale_val: int = None
    oreoSmallScale_des: str = None
    oreoTensScale_val: int = None
    oreoTensScale_des: str = None
    oreoHundsScale_val: int = None
    oreoHundsScale_des: str = None
    oreoLargeScale_val: int = None
    oreoLargeScale_des: str = None
    IUPAC_name: str = ''
    CAS_no: str = ''
 
    # Internal Processing Values
    mol: any = None
    molPixmap: any = None
    mwStr: str = None
    obStr: str = None
    isStr: str = None
    epStr: str = None
    eleList: any = None
    Catoms: int = None
    Hatoms: int = None
    Oatoms: int = None
    yoshidaMethod: str = 'Pfizer'

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
        self.molPixmap = pixmap
        return pixmap

    def mwFromMol(self):
        cmpdMW = Descriptors.MolWt(self.mol)
        self.MW = cmpdMW
        self.mwStr = "{:.2f}".format(cmpdMW)

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
        if self.EFG >= 1:
            self.oreo += 8
        elif self.EFG == 0:
            self.oreo += 1
        else:
            print("Error: O.R.E.O.S. EFG number gives unexpected value: " + EFG)   

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
        self.eleList = eleList
        return eleList

    def CHOFromEleComp(self):
        carbonAtoms = [ele[1] for ele in self.eleList if "C" in ele[0]]
        hydrogenAtoms = [ele[1] for ele in self.eleList if "H" in ele[0]]
        oxygenAtoms = [ele[1] for ele in self.eleList if "O" in ele[0]]
     
        # Get Number of Carbon Atoms
        if carbonAtoms == []:
            Catoms = 0
        elif carbonAtoms[0] == '':
            Catoms = 1
        else:
            Catoms = int(carbonAtoms[0])

        # Get Number of Hydrogen Atoms
        if hydrogenAtoms == []:
            Hatoms = 0
        elif hydrogenAtoms[0] == '':
            Hatoms = 1
        else:
            Hatoms = int(hydrogenAtoms[0])
        
        # Get Number of Oxygen Atoms
        if oxygenAtoms == []:
            Oatoms = 0
        elif oxygenAtoms[0] == '':
            Oatoms = 1
        else:
            Oatoms = int(oxygenAtoms[0])

        self.Catoms = Catoms
        self.Hatoms = Hatoms
        self.Oatoms = Oatoms
        return Catoms, Hatoms, Oatoms

    def RoSFromEleComp(self): 
        ruleSixCrit = 6*self.HEG - self.Catoms
        if ruleSixCrit > 0:
            ruleSix = " <b>(Explosive)</b>"
            self.oreo += 8
        else:
            ruleSix = " <b>(Not Explosive)</b>"
            self.oreo += 2 
        self.RoS_val = ruleSixCrit
        self.RoS_des = ruleSix

    def OBFromEleComp(self):
        oxygenBalance = (-1600*((2*self.Catoms)+(self.Hatoms/2)-self.Oatoms))/self.MW
        if oxygenBalance > 160:
            obRisk = "<b>(Low Risk)</b>"
            self.oreo += 2
        elif oxygenBalance > 80 and oxygenBalance <= 160:
            obRisk = "<b>(Medium Risk)</b>"
            self.oreo += 4
        elif oxygenBalance >= -120 and oxygenBalance <= 80:
            obRisk = "<b>(High Risk)</b>"
            self.oreo += 8
        elif oxygenBalance >= -240 and oxygenBalance < -120:
            obRisk = "<b>(Medium Risk)</b>"
            self.oreo += 4
        elif oxygenBalance < -240:
            obRisk = "<b>(Low Risk)</b>"
            self.oreo += 2
        self.OB_val = oxygenBalance
        self.OB_des = obRisk
        self.obStr = "{:.2f}".format(oxygenBalance)

    def oreoOnsetTadjustment(self):
        if self.onsetT != 'nan' and self.onsetT != '' and self.onsetT != None:
            print('Onset Temperature given as ' + str(self.onsetT))
            self.onsetT = float(self.onsetT)
            if self.onsetT <= 125:
                self.oreo += 8
            elif self.onsetT in range(126,200):
                self.oreo += 4
            elif self.onsetT in range(200,300):
                self.oreo += 2
            elif self.onsetT >= 300:
                self.oreo += 1
        else:
            self.onsetT = None
             
    def oreoSafeScaleCal(self):
        self.oreoSmallScale_val = self.oreo + 1
        self.oreoTensScale_val = self.oreo + 2
        self.oreoHundsScale_val = self.oreo + 4
        self.oreoLargeScale_val = self.oreo + 8

        scaleList = [self.oreoSmallScale_val, self.oreoTensScale_val, self.oreoHundsScale_val, self.oreoLargeScale_val]
        hazardList = []

        if self.onsetT == None:
            hazardValuesRangeList = [15, 22]
        elif self.onsetT != None:
            hazardValuesRangeList = [18, 28]

        for scale in scaleList:
            if scale < hazardValuesRangeList[0]:
                oreosHazard = "Low Hazard"
                hazardList.append(oreosHazard)
            elif scale in range(hazardValuesRangeList[0], hazardValuesRangeList[1]):  
                # Note: Remember the last integer in range isnt included thats why I can get away with a list of 2 items. 
                oreosHazard = "Medium Hazard"
                hazardList.append(oreosHazard)
            elif scale >= hazardValuesRangeList[1]:
                oreosHazard = "High Hazard"
                hazardList.append(oreosHazard)
            else:
                print('Error: ' + str(scale))
        print(hazardList)
        self.oreoSmallScale_des,self.oreoTensScale_des,self.oreoHundsScale_des,self.oreoLargeScale_des = hazardList

    def Td24FromThermalProps(self):
        # Estimation of Maximum Recommended Process Temperature To Avoid Hazardous Thermal Decomposition
        if self.initT != 'nan' and self.initT != '' and self.initT != None:
            print('Tinit given as ' + str(self.initT))
            self.initT = float(self.initT)
            d24Temp = 0.7*self.initT - 46
            print('T_D24 =' + str(d24Temp))

            self.Td24 = d24Temp

    def yoshidaFromThermalProps(self):
        # Determin if Original Yoshida calculation is to be used or if Pfizer modification is to be used.
        if self.yoshidaMethod == 'Yoshida':
            T_yoshida = self.onsetT
            ISConstant = 0.72
            EPConstant = 0.38
        elif self.yoshidaMethod == 'Pfizer':
            T_yoshida = self.initT
            ISConstant = 0.54
            EPConstant = 0.285

        # Perform the calculation
        if T_yoshida != 'nan' and T_yoshida != '' and T_yoshida != None and self.Q_dsc != 'nan' and self.Q_dsc != '' and self.Q_dsc != None:
                # Convert Qdsc to cal/g if needed
                print('Qdsc given as ' + str(self.Q_dsc) + ' ' + self.Qunits)
                self.Q_dsc = float(self.Q_dsc)
                if self.Qunits == 'J g<sup>-1</sup>' or self.Qunits == 'J g⁻¹':
                    calQ = self.Q_dsc/4.184
                elif self.Qunits == 'cal g<sup>-1</sup>' or self.Qunits == 'cal g⁻¹':
                    calQ = self.Q_dsc

                print(calQ)
             
                # Impact Sensitvity (IS)
                impactSens = log10(calQ) - ISConstant*(log10(abs(self.onsetT-25))) - 0.98
                print('IS = ' + str(impactSens))
             
                # Explosive Propagation (EP)
                exProp = log10(calQ) - EPConstant*(log10(abs(self.onsetT-25))) - 1.67
                print('EP = ' + str(exProp))
          
                # Warning Text
                if impactSens >= 0:
                    impact = ' <b>(Impact Sensitive)</b>'
                elif impactSens < 0:
                    impact = ' <b>(Not Impact Sensitive)</b>'
                if exProp >= 0:
                    explos = ' <b>(Propagates)</b>'
                elif exProp < 0:
                    explos = ' <b>(Should Not Propagate)</b>'

                # Save to thermalDexMolecule
                self.IS_val = impactSens
                self.IS_des = impact
                self.EP_val = exProp
                self.EP_des = explos
                self.isStr = "{:.2f}".format(impactSens)
                self.epStr = "{:.2f}".format(exProp)


    def IUPACNameFromSMILES(self):
        foundNames = get_compounds(self.SMILES, namespace='smiles')
        IUPACName = foundNames[0].iupac_name
        print(IUPACName)
        self.IUPAC_name = IUPACName

    def genCoreValues(self):
        self.genMol()
        self.mwFromMol()
        self.HEGFromMol()
        self.EFGFromMol()
        self.eleCompFromMol()
        self.CHOFromEleComp()
        self.RoSFromEleComp()
        self.OBFromEleComp()
        self.oreoOnsetTadjustment()
        self.oreoSafeScaleCal()
        self.yoshidaFromThermalProps()
        self.Td24FromThermalProps()

    def genAdditionalValues(self):
        self.IUPACNameFromSMILES()

    def genAllValues(self):
        self.genCoreValues()
        self.genAdditionalValues()

    def genCoreValuesFromMol(self):
        try:
            self.mwFromMol()
        except:
            window.showErrorMessage("Calculating MW from RDChem Mol.")

        try:
           self.HEGFromMol()
        except:
            window.showErrorMessage("Determining High Energy Groups from RDChem Mol.")

        try:
           self.EFGFromMol()
        except:
            window.showErrorMessage("Determining Explosive Groups from RDChem Mol.")

        try:
           self.eleCompFromMol()
        except:
            window.showErrorMessage("Determining Elemental Compostion from RDChem Mol.")
        
        try:
           self.CHOFromEleComp()
        except:
            window.showErrorMessage("Determining CHO from Elemental Compostion.")

        try:        
           self.RoSFromEleComp()
        except:
            window.showErrorMessage("Calculating Rule of Six from Elemental Compostion.")

        try:
           self.OBFromEleComp()
        except:
            window.showErrorMessage("Calculating Oxygen Balance from Elemental Compostion.")

        try:
           self.oreoOnsetTadjustment()
        except:
            window.showErrorMessage("Adjusting O.R.E.O.S. Calculation for Onset Temperature")

        try:
           self.oreoSafeScaleCal()
        except:
            window.showErrorMessage("Determining O.R.E.O.S. Hazard by Scale.")

        try:
           self.Td24FromThermalProps()
        except:
            window.showErrorMessage("Calculating T<sub>D24</sub> from Thermal Properties.")

        try:
           self.yoshidaFromThermalProps()
        except:
            window.showErrorMessage("Calculating Yoshida values from Thermal Properties.")





class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)

class QVLine(QFrame):
    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)

class ClickableLineEdit(QLineEdit):
    clicked = pyqtSignal() # signal when the text entry is left clicked

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton: self.clicked.emit()
        else: super().mousePressEvent(event)

class MolDrawer(QWidget):
    def __init__(self):
        super().__init__()

        self.init_ui()

    def init_ui(self):
        self.df = None
        #current_index = 0
        self.current_index = 0
        self.result_smiles = None
        self.error_flag = None
        # Set up the main layout
        layout = QVBoxLayout()

        # Create a tab widget
        tab_widget = QTabWidget()

        # Tab for molecule rendering
        molecule_tab = QWidget()
        molecule_layout = QVBoxLayout()

        # Display area for the molecular drawing
        self.mol_display = QGraphicsView(self)
        molecule_layout.addWidget(QLabel('Molecule:'))
        molecule_layout.addWidget(self.mol_display)

        ResultsContainLayout = QHBoxLayout()
        ResultsLeftLayout = QVBoxLayout()
        ResultsRightLayout = QVBoxLayout()

	#Add labels for calculated values
        self.mwText = 'MW: '
        self.HEGText = 'Number of High Energy Groups:'
        self.EFGText = 'Number of Explosive Functional Groups:'
        self.eleText = 'Elemental Composition: '
        self.RoSText = 'Rule of Six: '
        self.obText = 'Oxygen Balance: '
        self.ISText = 'Yoshida Impact Sensitivity: '
        self.EPText = 'Yoshida Explosive Propagation: '
        self.Td24Text = 'T<sub>D24</sub>: '

        self.mwLabel = QLabel(self.mwText)
        self.mwLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsLeftLayout.addWidget(self.mwLabel)
        self.HEGlabel = QLabel(self.HEGText)
        self.HEGlabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsLeftLayout.addWidget(self.HEGlabel)
        self.EFGlabel = QLabel(self.EFGText)
        self.EFGlabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsLeftLayout.addWidget(self.EFGlabel)
        self.eleLabel = QLabel(self.eleText)
        self.eleLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsLeftLayout.addWidget(self.eleLabel)
        self.RoSLabel = QLabel(self.RoSText)
        self.RoSLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.RoSLabel)
        self.obLabel = QLabel(self.obText)
        self.obLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.obLabel)
        self.ISLabel = QLabel(self.ISText)
        self.ISLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.ISLabel)
        self.EPLabel = QLabel(self.EPText)
        self.EPLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.EPLabel)
        self.Td24Label = QLabel(self.Td24Text)
        self.Td24Label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.Td24Label)

        ResultsContainLayout.addWidget(QVLine())
        ResultsContainLayout.addLayout(ResultsLeftLayout)
        ResultsContainLayout.addWidget(QVLine())
        ResultsContainLayout.addLayout(ResultsRightLayout)
        ResultsContainLayout.addWidget(QVLine())
        molecule_layout.addLayout(ResultsContainLayout)
        molecule_layout.addWidget(QHLine())

        self.tableLabel = QLabel('O.R.E.O.S. Assessment of Hazard by Scale:')
        molecule_layout.addWidget(self.tableLabel)
        tableLayout = QHBoxLayout()
        self.table = QTableWidget(1, 4)
        self.table.setHorizontalHeaderLabels(['<5 g', '5 to <100 g', '100 to 500 g', '>500 g'])
        self.table.verticalHeader().setVisible(False)
        #self.table.setStyleSheet("QTableWidget::item { border-bottom: 2px solid black; }")
        self.table.setMaximumHeight(53)
        self.table.setMaximumWidth(402)
        self.table.setMinimumHeight(53)
        self.table.setMinimumWidth(402)
        #self.table.setAlignment(Qt.AlignVCenter)
        tableLayout.addWidget(self.table)

        molecule_layout.addLayout(tableLayout)
        molecule_layout.addWidget(QHLine())

        # Input field for SMILES string
        #self.smiles_input = QLineEdit(self)
        self.smiles_input = ClickableLineEdit(self)
        self.smiles_input.clicked.connect(self.mrvToSMILES)
        molecule_layout.addWidget(QLabel('Enter SMILES String:'))
        molecule_layout.addWidget(self.smiles_input)

        InputContainLayout = QHBoxLayout()
        InputLeftLayout = QVBoxLayout()
        InputRightLayout = QVBoxLayout()
        numValidator = QRegExpValidator(QRegExp(r'[-]?\d+[.]?\d*'))
        posNumValidator = QRegExpValidator(QRegExp(r'\d+[.]?\d*'))

        # Input field for Name string
        self.name_input = QLineEdit(self)
        InputLeftLayout.addWidget(QLabel('Name:'))
        nameUnitsSubLayout = QHBoxLayout()
        nameUnitsSubLayout.addWidget(self.name_input)
        nameUnitLabel = QLabel('°C    ')
        nameUnitLabel.setStyleSheet('color: white')
        nameUnitsSubLayout.addWidget(nameUnitLabel)
        InputLeftLayout.addLayout(nameUnitsSubLayout)

        # Input field for mp string
        self.mp_input = QLineEdit(self)
        self.mp_input.setMaximumWidth(110)
        self.mp_input.setValidator(numValidator)
        self.mpEnd_input = QLineEdit(self)
        self.mpEnd_input.setMaximumWidth(110)
        self.mpEnd_input.setValidator(numValidator)
        InputLeftLayout.addWidget(QLabel('m.p.:'))
        mpUnitsSubLayout = QHBoxLayout()
        mpUnitsSubLayout.addWidget(self.mp_input)
        mpUnitsSubLayout.addWidget(QLabel('  to  '))
        mpUnitsSubLayout.addWidget(self.mpEnd_input)
        mpUnitsSubLayout.addWidget(QLabel('°C    '))
        InputLeftLayout.addLayout(mpUnitsSubLayout)

        # Input field for Q_DSC string
        self.Qdsc_input = QLineEdit(self)
        self.Qdsc_input.setValidator(posNumValidator)
        InputRightLayout.addWidget(QLabel('Q<sub>DSC</sub>:'))
        QUnitsSubLayout = QHBoxLayout()
        QUnitsSubLayout.addWidget(self.Qdsc_input)
        #QUnitsSubLayout.addWidget(QLabel('cal g<sup>-1</sup> '))
        self.QUnitsSelection = QComboBox(self)
        self.QUnitsSelection.addItems(['J g⁻¹', 'cal g⁻¹'])
        QUnitsSubLayout.addWidget(self.QUnitsSelection)
        InputRightLayout.addLayout(QUnitsSubLayout)

        # Input field for Onset string
        self.TE_input = QLineEdit(self)
        self.TE_input.setValidator(numValidator)
        InputRightLayout.addWidget(QLabel('Onset Temperature:'))
        TEUnitsSubLayout = QHBoxLayout()
        TEUnitsSubLayout.addWidget(self.TE_input)
        TEUnitsSubLayout.addWidget(QLabel('°C       '))
        InputRightLayout.addLayout(TEUnitsSubLayout)

        # Input field for Init string
        self.Tinit_input = QLineEdit(self)
        self.Tinit_input.setValidator(numValidator)
        InputRightLayout.addWidget(QLabel('Initiation Temperature:'))
        TinitUnitsSubLayout = QHBoxLayout()
        TinitUnitsSubLayout.addWidget(self.Tinit_input)
        TinitUnitsSubLayout.addWidget(QLabel('°C       '))
        InputRightLayout.addLayout(TinitUnitsSubLayout)

        # Input field for proj string
        self.proj_input = QLineEdit(self)
        InputLeftLayout.addWidget(QLabel('Project:'))
        projUnitsSubLayout = QHBoxLayout()
        projUnitsSubLayout.addWidget(self.proj_input)
        projUnitLabel = QLabel('°C    ')
        projUnitLabel.setStyleSheet('color: white')
        projUnitsSubLayout.addWidget(projUnitLabel)
        InputLeftLayout.addLayout(projUnitsSubLayout)

        InputContainLayout.addLayout(InputLeftLayout)
        #ResultsContainLayout.addWidget(QVLine())
        InputContainLayout.addLayout(InputRightLayout)
        #ResultsContainLayout.addWidget(QVLine())
        molecule_layout.addLayout(InputContainLayout)

        # Buttons
        buttonSpacerLabel = QLabel('hidden spacer')
        buttonSpacerLabel.setStyleSheet('color: white')
        molecule_layout.addWidget(buttonSpacerLabel)
        buttonContainLayout = QHBoxLayout()
        render_button = QPushButton('Evaluate Molecule', self)
        render_button.clicked.connect(self.render_molecule)
        render_button.setMaximumWidth(180)
        buttonContainLayout.addWidget(render_button)
        clear_button = QPushButton('Clear Sheet', self)
        clear_button.clicked.connect(self.resetToDefaultState)
        clear_button.setMaximumWidth(180)
        buttonContainLayout.addWidget(clear_button)
        #molecule_layout.addLayout(buttonContainLayout)
        msg_button = QPushButton('Error Message Test', self)
        msg_button.clicked.connect(self.errorTestingTool)
        msg_button.setMaximumWidth(180)
        buttonContainLayout.addWidget(msg_button)
        molecule_layout.addLayout(buttonContainLayout)

        molecule_tab.setLayout(molecule_layout)
        tab_widget.addTab(molecule_tab, "Add")


        # Tab for Search
        search_tab = QWidget()
        search_layout = QVBoxLayout()

        # Entry widget for searching
        lbl_search = QLabel('Search:')
        entry_search = QLineEdit()
        result_label = QLabel('click search')
        counter_label = QLabel('none')

        def search_database():
        #    query = entry_search.text().lower()
        #    results = [entry for entry in self.database if query in entry['Name'].lower() or query in entry['Formula'].lower()]

        #    list_search_results.clear()
        #    for result in results:
        #        list_search_results.addItem(f"{result['Name']} ({result['Formula']})")

             self.df = pd.read_csv(defaultDB, encoding='mbcs')
             show_result(self)


        def show_result(self):
             #print(self)
             layout = self.layout()
             if self.error_flag is not None:
                  self.error_message.setText('')
                  layout.removeWidget(self.error_message)
                  self.error_flag = None
             if self.df is not None and not self.df.empty:
                  #print(self.current_index)
                  current_row = self.df.iloc[self.current_index]
                  #print(self.result_smiles)
                  self.result_smiles = current_row['Molecule']
                  #print(self.result_smiles)
                  result_text = f"SMILES: {current_row['Molecule']}\nName: {current_row['Name']}\nHEG: {current_row['Number of High Energy Groups']}\nmp: {current_row['Melting point']}\nMW: {current_row['Molecular Weight']}\nTE: {current_row['Thermal Event']}\nProject: {current_row['Project']}"
                  result_label.setText(result_text)
                  result_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                  result_label.setWordWrap(True) 
                  counter_label.setText(f"Result {self.current_index + 1} of {len(self.df)}")
                  search_layout.addWidget(prev_button)
                  search_layout.addWidget(next_button)
                  #print(self.current_index)

                  if self.result_smiles is not None:
                      mol = MolFromSmiles(self.result_smiles)

                  else: 
                      mol = None

                  if mol is not None:
                      # Generate a molecular drawing as a PNG image
                      img = Draw.MolToImage(mol)

                      # Convert the image to a byte array
                      byte_array = BytesIO()
                      img.save(byte_array, format='PNG')

                      # Convert the byte array to a QPixmap and display it
                      pixmap = QPixmap()
                      pixmap.loadFromData(byte_array.getvalue())
                      scene = QGraphicsScene()
                      scene.addPixmap(pixmap)
                      self.mol_result_display.setScene(scene)

        def prev_result(self):
             if self.current_index > 0:
                  #search_layout.removeWidget(self.mol_result_display)
                  #search_layout.removeWidget(self.molLabel)
                  self.current_index -= 1
                  show_result(self)

        def next_result(self):
             if self.df is not None:
                  if self.current_index < len(self.df) - 1:
                       #search_layout.removeWidget(self.mol_result_display)
                       #search_layout.removeWidget(self.molLabel)
                       self.current_index += 1
                       show_result(self)


        # Search Buttons & display area for the molecular drawing
        self.mol_result_display = QGraphicsView(self)
        self.molLabel = QLabel('Molecule:')
        search_layout.addWidget(self.molLabel)
        search_layout.addWidget(self.mol_result_display)

        btn_search = QPushButton('Search', clicked=search_database)
        prev_button = QPushButton('Previous')
        next_button = QPushButton('Next')
        prev_button.clicked.connect(lambda: prev_result(self))
        next_button.clicked.connect(lambda: next_result(self))

        search_layout.addWidget(lbl_search)
        search_layout.addWidget(entry_search)
        search_layout.addWidget(result_label)
        search_layout.addWidget(counter_label)
        search_layout.addWidget(btn_search)

        search_tab.setLayout(search_layout)
        tab_widget.addTab(search_tab, "Search")


        # Tab for "About" information

        #def hover(url):
        #    if url:
        #        QToolTip.showText(QCursor.pos(), titles.get(url, url))
        #    else:
        #        QToolTip.hideText()

        about_tab = QWidget()
        about_layout = QVBoxLayout()
        about_title = QLabel("<b>About ThermalDex</b>\n\n")
        about_blank = QLabel("\nVersion: " + versionNumber + " (This is currently an alpha build)\n")
        about_text = QLabel("\n\nThis is a simple tool for assessing and recording the potential thermal hazards assoicated with a molecule. It uses the <b>'O.R.E.O.S.'</b> assement scale and other ideas that can be read about in <a href=\"https://pubs.acs.org/doi/10.1021/acs.oprd.0c00467\"><em>Org. Process Res. Dev.</em> 2021, 25, 2, 212–224</a> by Jeffrey B. Sperry et. al.")
        iconLabel = QLabel()
        iconImage = QPixmap(".\\_core\\ThermalDexIcon.jpg")
        scaledIcon = iconImage.scaled(400, 400, Qt.KeepAspectRatio)
        iconLabel.setText("test") #.setPixmap(scaledIcon)
        
        scene = QGraphicsScene()
        pic = QGraphicsPixmapItem()
        pic.setPixmap(scaledIcon) #QPixmap.fromImage(iconImage))
        scene.setSceneRect(0, 0, 400, 400)
        scene.addItem(pic)
        view = QGraphicsView()
        view.setScene(scene)
        view.setStyleSheet("border: 0px")
        altNames_text = QLabel("The alternative name for this project is <b>C.O.O.K.I.E.S.</b> which stands for <em>C</em>alculation <em>O</em>f <em>O</em>.R.E.O.S., <em>K</em>illing <em>I</em>nelegant <em>E</em>xcel <em>S</em>olutions. This would be the offical relase name, but it felt a bit glib.")
        moreabout_text = QLabel("This tool has been developed by Tom Sheridan and Matt Mulheir using PyQt5 and RDKit. You can contribute to this project on <a href=\"https://github.com/Tomisall/ProcessSafteyDB\">GitHub</a>.")
        #about_text.setAlignment(Qt.AlignCenter)
        about_layout.addWidget(about_title)
        about_layout.addWidget(about_text)
        about_layout.addWidget(about_blank)
        about_layout.addWidget(view) #iconLabel)
        about_layout.addWidget(altNames_text)
        about_layout.addWidget(moreabout_text)
        about_text.setTextFormat(Qt.RichText)
        about_text.setOpenExternalLinks(True)
        about_text.setTextInteractionFlags(Qt.LinksAccessibleByMouse)
        #about_text.linkHovered.connect(hover)
        about_text.setWordWrap(True) 
        altNames_text.setWordWrap(True) 
        moreabout_text.setTextInteractionFlags(Qt.LinksAccessibleByMouse)
        moreabout_text.setOpenExternalLinks(True)
        moreabout_text.setWordWrap(True)

        #abouttestLabel = QLabel("<a href=\"http://www.google.com\">'Click this link to go to Google'</a>")
        #about_layout.addWidget(abouttestLabel)
        #abouttestLabel.setOpenExternalLinks(True)
        #abouttestLabel.setTextInteractionFlags(Qt.LinksAccessibleByMouse)

        about_tab.setLayout(about_layout)
        tab_widget.addTab(about_tab, "About")

        layout.addWidget(tab_widget)



        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 600, 825)
        self.setWindowTitle('ThermalDex')

    def mrvToSMILES(self):
        try:
            copyMolXML = pyperclip.paste()
            copyMolReal = rdmolfiles.MolFromMrvBlock(copyMolXML)
            smilesFromMarvin = rdmolfiles.MolToSmiles(copyMolReal)
            print(smilesFromMarvin)
            pyperclip.copy(smilesFromMarvin)

        except:
            print("No mrv XML found.")

    def showErrorMessage(self, errorCode):
        self.msg = QMessageBox()
        self.msg.setIcon(QMessageBox.Warning)
        self.msg.setText("An error has occured")
        self.msg.setInformativeText("The source of the problem seems to be with: " + errorCode + "\nTry again and contact developer if the problem persists.")
        self.msg.setWindowTitle("ThermalDex Error")
        self.msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        returnValue = self.msg.exec()

    def errorTestingTool(self):
        self.showErrorMessage("This is an Alternative test method for error handling")

    def clearTheCalcdValues(self):
        scene = QGraphicsScene()
        self.mol_display.setScene(scene)
        self.mwLabel.setText(self.mwText)
        self.HEGlabel.setText(self.HEGText)
        self.EFGlabel.setText(self.EFGText)
        self.eleLabel.setText(self.eleText)
        self.RoSLabel.setText(self.RoSText)
        self.obLabel.setText(self.obText)
        self.ISLabel.setText(self.ISText)
        self.EPLabel.setText(self.EPText)
        self.Td24Label.setText(self.Td24Text)
        clearEntry = ['', '', '', '']
        for i, entry in enumerate(clearEntry):
            clear = QTableWidgetItem(entry)
            self.table.setItem(0, i, clear)
            clear.setBackground(QColor(255, 255, 255))

    def clearUserValues(self):
        self.smiles_input.setText('')
        self.name_input.setText('')
        self.mp_input.setText('')
        self.mpEnd_input.setText('')
        self.Qdsc_input.setText('')
        self.TE_input.setText('')
        self.Tinit_input.setText('')
        self.proj_input.setText('')

    def resetToDefaultState(self):
        self.clearTheCalcdValues()
        self.clearUserValues()
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)

    def getColorForValue(self, hazardClass):
        # Example color-coding logic
        if hazardClass == 'High Hazard':
            return QColor(255, 0, 0)  # Red
        elif hazardClass == 'Medium Hazard':
            return QColor(255, 255, 0)  # Yellow
        elif hazardClass == 'Low Hazard':
            return QColor(0, 255, 0)  # Green
        else:
            return QColor(0, 0, 255)  # Blue

    def render_molecule(self):
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)

        # Get the SMILES string from the input field
        smiles = self.smiles_input.text()
        name = self.name_input.text()
        mp = self.mp_input.text()
        Qdsc = self.Qdsc_input.text()
        QUnits = self.QUnitsSelection.currentText()
        TE = self.TE_input.text()
        Tinit = self.Tinit_input.text()
        proj = self.proj_input.text()

        # To be removed and replaced by call to thermalDexMolecule, see later.
        writeSmiles = '"'+smiles+'"' #repr(str(smiles))
        writeName = '"'+name+'"'
        writemp = '"'+mp+'"'
        writeTE = '"'+TE+'"'
        writeProj = '"'+proj+'"'

        # Create an RDKit molecule from the SMILES string
        addedMolecule = thermalDexMolecule(SMILES=smiles, name=name, mp=mp, Q_dsc=Qdsc, Qunits=QUnits, onsetT=TE, initT=Tinit, proj=proj)
        if smiles != '' and smiles is not None:
            addedMolecule.mol = MolFromSmiles(smiles)

        else:
            addedMolecule.mol = None
        
        if addedMolecule.mol is not None:
            # Make Pixmap Image to Display.
            pixmap = addedMolecule.molToQPixmap()
            scaledPixmap = pixmap.scaled(300, 300, Qt.KeepAspectRatio)
            scene = QGraphicsScene()
            #scene.setSceneRect(0, 0, 400, 400)
            scene.addPixmap(scaledPixmap) #pixmap)
            self.mol_display.setScene(scene)

            # Calculate Core Properties
            addedMolecule.genCoreValuesFromMol()

            # Format and Display Properties
            self.mwLabel.setText('MW: ' + addedMolecule.mwStr)
            self.HEGlabel.setText('Number of High Energy Groups: ' + str(addedMolecule.HEG))
            self.EFGlabel.setText('Number of Explosive Functional Groups: ' + str(addedMolecule.EFG)) 
            self.eleLabel.setText('Elemental Composition: ' + addedMolecule.eleComp)
            self.RoSLabel.setText('Rule of Six: ' + str(addedMolecule.RoS_val) + addedMolecule.RoS_des)
            self.obLabel.setText('Oxygen Balance: ' + addedMolecule.obStr + ' ' + addedMolecule.OB_des)
            self.table.clearContents()
            hazardList = [addedMolecule.oreoSmallScale_des, addedMolecule.oreoTensScale_des, addedMolecule.oreoHundsScale_des, addedMolecule.oreoLargeScale_des]
            # Add values to cells
            for i, hazardClass in enumerate(hazardList):
                item = QTableWidgetItem(hazardClass)
                self.table.setItem(0, i, item)

                # Color code cells based on values
                classColor = self.getColorForValue(hazardClass)
                print(classColor)
                item.setBackground(classColor)

            if addedMolecule.isStr != None:
                self.ISLabel.setText('Yoshida Impact Sensitivity: ' + addedMolecule.isStr + addedMolecule.IS_des)
            if addedMolecule.epStr != None:
                self.EPLabel.setText('Yoshida Explosive Propagation: ' + addedMolecule.epStr + addedMolecule.EP_des)
            if addedMolecule.Td24 != '' and addedMolecule.Td24 != 'nan' and addedMolecule.Td24 != None:
                d24Str = "{:.1f}".format(addedMolecule.Td24)
                self.Td24Label.setText('T<sub>D24</sub>: ' + '<b>' + d24Str + ' °' + 'C' + '</b>')
            if addedMolecule.onsetT != 'nan' and addedMolecule.onsetT != '' and addedMolecule.onsetT != None and addedMolecule.onsetT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Onset Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 201
            if addedMolecule.initT != 'nan' and addedMolecule.initT != '' and addedMolecule.initT != None and addedMolecule.initT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Initiation Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 202

            addMol = open(defaultDB, 'a')
            #addMol.write(writeSmiles + ',' + writeName + ',' + HEG + ',' + writemp + ',' + mwStr + ',' + writeTE + ',' + writeProj + '\n')
            print('\n')
            print(addedMolecule)
            print('\n')


        else:
            # Display an error message if the SMILES string is invalid
            self.error_message = QLabel('Invalid SMILES string. Please enter a valid SMILES.')
            #layout = self.layout()
            layout.addWidget(self.error_message)
            self.error_flag = 100

if __name__ == '__main__':
    defaultDB, highEnergyGroups, expEnergyGroups = importConfig()
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('.\\_core\\ThermalDexIcon.ico'))
    window = MolDrawer()
    window.show()
    sys.exit(app.exec_())
