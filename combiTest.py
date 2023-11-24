import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTabWidget
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Mol
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt
from io import BytesIO
import pandas as pd
import re

def importConfig():
    conf = open('ThermalDex.config', 'r')
    confCounter = 0
    for line in conf:
        #print(confCounter)
        if confCounter == 4:
           defaultDB = line.strip("\n")
           confCounter += 1
        elif confCounter == 8:
           highEnergyGroups = line.strip("\n")
           confCounter += 1
        else:
           confCounter += 1

    #print(defaultDB)
    return defaultDB, highEnergyGroups

class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)

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

	#Add labels for calculated values
        self.mwLabel = QLabel('MW: ')
        self.mwLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        molecule_layout.addWidget(self.mwLabel)
        self.HEGlabel = QLabel('Number of High Energy Groups:')
        self.HEGlabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        molecule_layout.addWidget(self.HEGlabel)

        molecule_layout.addWidget(QHLine())

        # Input field for SMILES string
        self.smiles_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Enter SMILES String:'))
        molecule_layout.addWidget(self.smiles_input)

        # Input field for Name string
        self.name_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Name:'))
        molecule_layout.addWidget(self.name_input)

        # Input field for mp string
        self.mp_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('m.p.:'))
        molecule_layout.addWidget(self.mp_input)

        # Input field for TE string
        self.TE_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Thermal Event:'))
        molecule_layout.addWidget(self.TE_input)

        # Input field for proj string
        self.proj_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Project:'))
        molecule_layout.addWidget(self.proj_input)

        # Button to render the molecule
        render_button = QPushButton('Render Molecule', self)
        render_button.clicked.connect(self.render_molecule)
        molecule_layout.addWidget(render_button)

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
                  counter_label.setText(f"Result {self.current_index + 1} of {len(self.df)}")
                  search_layout.addWidget(prev_button)
                  search_layout.addWidget(next_button)
                  #print(self.current_index)

                  if self.result_smiles is not None:
                      mol = Chem.MolFromSmiles(self.result_smiles)

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
        about_tab = QWidget()
        about_layout = QVBoxLayout()

        about_text = QLabel("This is a simple molecule drawer using PyQt5 and RDKit.")
        #about_text.setAlignment(QtCore.Qt.AlignCenter)
        about_layout.addWidget(about_text)

        about_tab.setLayout(about_layout)
        tab_widget.addTab(about_tab, "About")

        layout.addWidget(tab_widget)



        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 400, 700)
        self.setWindowTitle('ThermalDex')

    def render_molecule(self):
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)
        # Get the SMILES string from the input field
        smiles = self.smiles_input.text()
        name = self.name_input.text()
        mp = self.mp_input.text()
        TE = self.TE_input.text()
        proj = self.proj_input.text()
        writeSmiles = '"'+smiles+'"' #repr(str(smiles))
        writeName = '"'+name+'"'
        writemp = '"'+mp+'"'
        writeTE = '"'+TE+'"'
        writeProj = '"'+proj+'"'

        # Create an RDKit molecule from the SMILES string
        if smiles != '' and smiles is not None:
            mol = Chem.MolFromSmiles(smiles)

        else:
            mol = None
        
        test = 'hi'

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
            self.mol_display.setScene(scene)
            cmpdMW = Descriptors.MolWt(mol)
            mwStr = "{:.2f}".format(cmpdMW)
            self.mwLabel.setText('MW: ' + mwStr)
            fullMatch = 0
            with open("HighEnergyGroups.csv", "r") as HEGroups:
               for line in HEGroups:
                    #print(line)
                    HeSubstructure = Chem.MolFromSmiles(line)
                    fullmatchList = Mol.GetSubstructMatches(mol, HeSubstructure)
                    fullMatch += len(fullmatchList)

            HEG = str(fullMatch)
            self.HEGlabel.setText('Number of High Energy Groups: ' + HEG)
            chemForm = rdMolDescriptors.CalcMolFormula(mol)
            print(chemForm)
            match = re.findall(r"[A-Z][a-z]?\d*|\((?:[^()]*(?:\(.*\))?[^()]*)+\)\d+", chemForm, re.I) #re.findall(r"([a-z]+)([0-9]+)", chemForm, re.I) #match
            if match:
               items = match #match.groups()
            print(items)
            addMol = open(defaultDB, 'a')
            addMol.write(writeSmiles + ',' + writeName + ',' + HEG + ',' + writemp + ',' + mwStr + ',' + writeTE + ',' + writeProj + '\n')


        else:
            # Display an error message if the SMILES string is invalid
            self.error_message = QLabel('Invalid SMILES string. Please enter a valid SMILES.')
            #layout = self.layout()
            layout.addWidget(self.error_message)
            self.error_flag = 100

if __name__ == '__main__':
    defaultDB, highEnergyGroups = importConfig()
    app = QApplication(sys.argv)
    window = MolDrawer()
    window.show()
    sys.exit(app.exec_())
