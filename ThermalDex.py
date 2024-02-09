#import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTableWidget, QTableWidgetItem, QTabWidget, QGraphicsPixmapItem,  QMessageBox, QComboBox #, QTableView, QToolTip
from PyQt5.QtGui import QPixmap, QColor, QIcon, QRegExpValidator #QValidator #, QCursor
from PyQt5.QtCore import Qt, QRegExp, pyqtSignal
#from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles
#from io import BytesIO
#from numpy import log10
#from pubchempy import get_compounds
#from dataclasses import dataclass, field, asdict
from thermDex.thermDexMolecule import * #thermalDexMolecule
from thermDex.thermDexReport import *
from thermDex.thermDexHTMLRep import *
from thermDex.attachedFileManager import *
#clea
#import re
import pyperclip
from pandasTests import *
import configparser

versionNumber = "0.8.6"

try:
    import pyi_splash
    pyi_splash.close()
    window.raise_()
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

def altImportConfig():
    config = configparser.ConfigParser()
    config.read('./_core/ThermalDex.ini')
    defaultDB = config.get('Database', 'defaultDB')
    highEnergyGroups = config.get('Lists of Groups', 'highEnergyGroups')
    expEnergyGroups = config.get('Lists of Groups', 'expEnergyGroups')
    yosidaMethod = config.get('Calculations', 'yoshidaCalcs')
    return defaultDB, highEnergyGroups, expEnergyGroups, yosidaMethod

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
        self.selectedDatabase = None
        self.error_flag = None
        # Set up the main layout
        layout = QVBoxLayout()

        # Create a tab widget
        self.tab_widget = QTabWidget()

        # Tab for molecule rendering
        self.molecule_tab = QWidget()
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

        # Input field for Hammer Drop Test
        hamSubLayout = QHBoxLayout()
        InputLeftLayout.addWidget(QLabel('Hammer Drop Test:'))
        self.hamSelection = QComboBox(self)
        self.hamSelection.addItems(['Not Performed', 'Positive', 'Negative'])
        hamSubLayout.addWidget(self.hamSelection)
        hamUnitLabel = QLabel('°C    ')
        hamUnitLabel.setStyleSheet('color: white')
        hamSubLayout.addWidget(hamUnitLabel)
        InputLeftLayout.addLayout(hamSubLayout)

        # Input field for Friction Test
        fricSubLayout = QHBoxLayout()
        InputRightLayout.addWidget(QLabel('Friction Test:'))
        self.fricSelection = QComboBox(self)
        self.fricSelection.addItems(['Not Performed', 'Positive', 'Negative'])
        fricSubLayout.addWidget(self.fricSelection)
        fricUnitLabel = QLabel('°C    ')
        fricUnitLabel.setStyleSheet('color: white')
        fricSubLayout.addWidget(fricUnitLabel)
        InputRightLayout.addLayout(fricSubLayout)

        InputContainLayout.addLayout(InputLeftLayout)
        #ResultsContainLayout.addWidget(QVLine())
        InputContainLayout.addLayout(InputRightLayout)
        #ResultsContainLayout.addWidget(QVLine())
        molecule_layout.addLayout(InputContainLayout)

        # Attach Data Files
        filesSubLayout = QHBoxLayout()
        self.attachedFilesLabel = QLabel('Attached Files:')
        self.attachedFilesLabel.resize(120, 120)
        filesSubLayout.addWidget(self.attachedFilesLabel)
        self.filesCount = QLabel('0 Attached Files')
        self.filesCount.resize(90, 120)
        filesSubLayout.addWidget(self.filesCount)
        self.attach_button = QPushButton('Add/View Files', self)
        self.attach_button.clicked.connect(self.openFileManager) 
        self.attach_button.setMaximumWidth(140)
        filesSubLayout.addWidget(self.attach_button)
        attachSpacerLabel = QLabel('hidden spacer')
        attachSpacerLabel.setStyleSheet('color: white')
        filesSubLayout.addWidget(attachSpacerLabel)
        attachSpacerLabelTwo = QLabel('hidden spacer')
        attachSpacerLabelTwo.setStyleSheet('color: white')
        filesSubLayout.addWidget(attachSpacerLabelTwo)
        molecule_layout.addSpacing(15)
        molecule_layout.addLayout(filesSubLayout)

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
        msg_button = QPushButton('Save to PDF', self)
        msg_button.clicked.connect(self.createReport)
        msg_button.setMaximumWidth(180)
        buttonContainLayout.addWidget(msg_button)
        molecule_layout.addLayout(buttonContainLayout)

        self.molecule_tab.setLayout(molecule_layout)
        self.tab_widget.addTab(self.molecule_tab, "Add")

        # Hide File Handling Widgets 
        self.attachedFilesLabel.hide()
        self.filesCount.hide()
        self.attach_button.hide()

        # Tab for Search
        search_tab = QWidget()
        search_layout = QVBoxLayout() #molecule_layout #

        # Entry widget for searching
        lbl_search = QLabel('Search:')
        entry_search = QLineEdit()
        searchTypeSelection = QComboBox(self)
        searchTypeSelection.addItems(['SMILES (substructure)', 'Name', 'Project', 'Compounds with smaller Qdsc', 'Compounds with larger Qdsc', 'Compounds with smaller Tinit', 'Compounds with larger Tinit', 'Compounds with smaller Tonset', 'Compounds with larger Tonset', 'Compounds with smaller Td24', 'Compounds with larger Td24', 'Compounds with smaller O.R.E.O.S. score at >500 g', 'Compounds with larger O.R.E.O.S. score at >500 g'])
        result_label = QLabel('click search')
        counter_label = QLabel('none')
        view_test_button = QPushButton('edit', self)
        view_test_button.clicked.connect(self.changeTabForEditing) 
        search_layout.addWidget(view_test_button)


        def search_database():
             self.readDatabase = pd.read_csv(defaultDB) #, index_col=0) #, encoding='mbcs')
             searchTest = MolFromSmiles(entry_search.text())
             print(entry_search.text())

             if entry_search.text() != '' and entry_search.text() != None and searchTest is not None:
                 smilesList = self.readDatabase['SMILES'].tolist()  #.index.values
                 print(smilesList)
                 foundList = []
                 for smile in smilesList:
                     searchStructure = MolFromSmiles(smile)
                     fullmatchList = Mol.GetSubstructMatches(searchStructure, searchTest)
                     if len(fullmatchList) > 0:
                         print('Substructure Match Found: ' + smile)
                         foundList += [smile]
                 print('\n\n')
                 print(foundList)

                 indexList = []
                 for foundMatch in foundList:
                     row_index = self.readDatabase.index[self.readDatabase['SMILES'] == foundMatch].tolist()
                     indexList += row_index

                 print(indexList)
                 print('\n\n')
                 foundDataFrame = self.readDatabase.iloc[indexList]
                 print(foundDataFrame)
                 print('\n\n')

                 foundDataFrame.reset_index(drop=True)

                 self.selectedDatabase = foundDataFrame
                 show_result(self, foundDataFrame, True)

             elif entry_search.text() == '' or entry_search.text() == None:    
                 self.selectedDatabase = self.readDatabase
                 show_result(self, self.readDatabase, True)

             else:
                 otherSearchIndex = []
                 otherSearchIndex += [self.readDatabase.isin([entry_search.text()]).any(axis=1).idxmax()]
                 print(otherSearchIndex)
                 otherFoundDataFrame = self.readDatabase.iloc[otherSearchIndex]
                 self.selectedDatabase = otherFoundDataFrame
                 show_result(self, otherFoundDataFrame, True)
                 if otherFoundDataFrame.empty:
                      errorInfo = "No matches found"
                      interactiveErrorMessage(self, errorInfo)
      


        def show_result(self, Database, resetIndex):
             #print(self)
             layout = self.layout()
             if self.error_flag is not None:
                  self.error_message.setText('')
                  layout.removeWidget(self.error_message)
                  self.error_flag = None
             if resetIndex == True:
                  self.current_index = 0             
             if Database is not None and not Database.empty:
                  current_row = Database.iloc[self.current_index]
                  print(current_row)
                  print('\n\n')
                  dictRow = current_row.to_dict()
                  print(dictRow)
                  print('\n\n')
                  readMolecule = thermalDexMolecule(**dictRow)
                  print(readMolecule.MW)
                  self.result_smiles = current_row['SMILES']
                  result_text = f"SMILES: {current_row['SMILES']}\nName: {current_row['name']}\nHEG: {current_row['HEG']}\nmp: {current_row['mp']} to {current_row['mpEnd']}\nMW: {current_row['MW']}\nT<sub>onset</sub>: {current_row['onsetT']}\nT<sub>init</sub>: {current_row['initT']}\nProject: {current_row['proj']}"
                  result_label.setText(result_text)
                  result_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                  result_label.setWordWrap(True) 
                  counter_label.setText(f"Result {self.current_index + 1} of {len(Database)}")
                  search_layout.addWidget(prev_button)
                  search_layout.addWidget(next_button)

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
                  show_result(self, self.selectedDatabase, False)

        def next_result(self, Database):
             if Database is not None:
                  if self.current_index < len(Database) - 1:
                       #search_layout.removeWidget(self.mol_result_display)
                       #search_layout.removeWidget(self.molLabel)
                       self.current_index += 1
                       show_result(self, self.selectedDatabase, False)

        # Search Buttons & display area for the molecular drawing
        self.mol_result_display = QGraphicsView(self)
        self.molLabel = QLabel('Molecule:')
        search_layout.addWidget(self.molLabel)
        search_layout.addWidget(self.mol_result_display)

        btn_search = QPushButton('Search', clicked=search_database)
        prev_button = QPushButton('Previous')
        next_button = QPushButton('Next')
        prev_button.clicked.connect(lambda: prev_result(self))
        next_button.clicked.connect(lambda: next_result(self, self.selectedDatabase))

        search_layout.addWidget(lbl_search)
        search_layout.addWidget(entry_search)
        search_layout.addWidget(searchTypeSelection)
        search_layout.addWidget(result_label)
        search_layout.addWidget(counter_label)
        search_layout.addWidget(btn_search)

        search_tab.setLayout(search_layout)
        self.tab_widget.addTab(search_tab, "Search")


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
        self.tab_widget.addTab(about_tab, "About")

        layout.addWidget(self.tab_widget)



        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 600, 875)
        self.setWindowTitle('ThermalDex')

    def changeTabForEditing(self):
        try:
            editDB = self.selectedDatabase.fillna('')
            current_row = editDB.iloc[self.current_index]
            dictRow = current_row.to_dict()
            print(dictRow)
            print('\n\n')
            readMolecule = thermalDexMolecule(**dictRow)
            readMolecule.genMol()
            # Make Pixmap Image to Display.
            pixmap = readMolecule.molToQPixmap()
            scaledPixmap = pixmap #.scaled(550, 275, Qt.KeepAspectRatio)
            scene = QGraphicsScene()
            #scene.setSceneRect(0, 0, 400, 400)
            scene.addPixmap(scaledPixmap) #pixmap)
            self.mol_display.setScene(scene)
            readMolecule.molPixmap = None
            self.smiles_input.setText(readMolecule.SMILES)
            self.name_input.setText(readMolecule.name)
            self.mp_input.setText(str(readMolecule.mp))
            self.mpEnd_input.setText(str(readMolecule.mpEnd))
            self.Qdsc_input.setText(str(readMolecule.Q_dsc))
            comboIndex = self.QUnitsSelection.findText(readMolecule.Qunits)
            self.QUnitsSelection.setCurrentIndex(comboIndex)
            self.TE_input.setText(str(readMolecule.onsetT))
            self.Tinit_input.setText(str(readMolecule.initT))
            self.proj_input.setText(str(readMolecule.proj))
            niceMWStr = "{:.2f}".format(readMolecule.MW)
            self.mwLabel.setText('MW: ' + niceMWStr) #str(readMolecule.MW))
            self.HEGlabel.setText('Number of High Energy Groups: ' + str(readMolecule.HEG))
            self.EFGlabel.setText('Number of Explosive Functional Groups: ' + str(readMolecule.EFG)) 
            self.eleLabel.setText('Elemental Composition: ' + readMolecule.eleComp)
            self.RoSLabel.setText('Rule of Six: ' + str(readMolecule.RoS_val) + readMolecule.RoS_des)
            self.obLabel.setText('Oxygen Balance: ' + str(readMolecule.obStr) + ' ' + readMolecule.OB_des)
            self.table.clearContents()
            hazardList = [readMolecule.oreoSmallScale_des, readMolecule.oreoTensScale_des, readMolecule.oreoHundsScale_des, readMolecule.oreoLargeScale_des]
            # Add values to cells
            for i, hazardClass in enumerate(hazardList):
                item = QTableWidgetItem(hazardClass)
                self.table.setItem(0, i, item)
   
                # Color code cells based on values
                classColor = self.getColorForValue(hazardClass)
                print(classColor)
                item.setBackground(classColor)
            
            fileCounter = self.countFiles(defaultDB)
            self.filesCount.setText(f"{str(fileCounter)} Attached Files")
            self.attachedFilesLabel.show()
            self.filesCount.show()
            self.attach_button.show()
        
            try:
                isStr = "{:.2f}".format(readMolecule.IS_val)
            except:
                isStr = ''
            self.ISLabel.setText('Yoshida Impact Sensitivity: ' + isStr + readMolecule.IS_des)
            try:
                epStr = "{:.2f}".format(readMolecule.EP_val)
            except:
                epStr = ''
            self.EPLabel.setText('Yoshida Explosive Propagation: ' + epStr + readMolecule.EP_des)
            try:
                d24Str = "{:.1f}".format(readMolecule.Td24)
                self.Td24Label.setText('T<sub>D24</sub>: <b>' + d24Str + ' °C</b>')        
            except:
                pass
            self.tab_widget.setCurrentWidget(self.molecule_tab)   #0) #self.tab_widget.findChild(QWidget, "Add"))

        except:
            window.showErrorMessage("Editing current selection. Has a molecule been found?")

    def mrvToSMILES(self):
        try:
            copyMolXML = pyperclip.paste()
            copyMolReal = rdmolfiles.MolFromMrvBlock(copyMolXML)
            smilesFromMarvin = rdmolfiles.MolToSmiles(copyMolReal)
            print(smilesFromMarvin)
            pyperclip.copy(smilesFromMarvin)

        except:
            print("No mrv XML found.")

    def openFileManager(self):
        self.fileWindow = FileManagementWindow(self.smiles_input.text(), defaultDB, window)
        self.fileWindow.show()
        self.fileWindow.raise_()
        self.fileWindow.activateWindow()

    def findFolder(self, database):
        try:
            storedData = pd.read_csv(database)
            row_index = storedData.index[storedData['SMILES'] == self.smiles_input.text()].tolist()
            folderInfo = storedData['dataFolder'][row_index[0]]
            return folderInfo

        except:
            return ''
     

    def countFiles(self, database):
        folderInfo = self.findFolder(database)
        files = [f for f in listdir(folderInfo) if isfile(join(folderInfo, f))]
        fileCount = len(files)
        return fileCount

    def showErrorMessage(self, errorCode):
        self.msg = QMessageBox()
        self.msg.setIcon(QMessageBox.Warning)
        self.msg.setText("An error has occured")
        self.msg.setInformativeText("The source of the problem seems to be with: " + errorCode + "\nTry again and contact developer if the problem persists.")
        self.msg.setWindowTitle("ThermalDex Error")
        self.msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        returnValue = self.msg.exec()

    def interactiveErrorMessage(self, errorInfo):
        self.interactMsg = QMessageBox()
        self.interactMsg.setIcon(QMessageBox.Information)
        self.interactMsg.setText("Action Needed")
        self.interactMsg.setInformativeText(errorInfo)
        self.interactMsg.setWindowTitle("ThermalDex - Info Box")
        self.interactMsg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        returnValue = self.interactMsg.exec()
        return returnValue

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
        self.hamSelection.setCurrentIndex(0)
        self.fricSelection.setCurrentIndex(0)
        self.attachedFilesLabel.hide()
        self.filesCount.hide()
        self.attach_button.hide()

    def resetToDefaultState(self):
        self.clearTheCalcdValues()
        self.clearUserValues()
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)

    def genCoreValuesFromMol(self, molecule):
        try:
            molecule.mwFromMol()
        except:
            window.showErrorMessage("Calculating MW from RDChem Mol.")

        try:
            molecule.HEGFromMol(highEnergyGroups)
        except:
            window.showErrorMessage("Determining High Energy Groups from RDChem Mol.")

        try:
           molecule.EFGFromMol(expEnergyGroups)
        except:
            window.showErrorMessage("Determining Explosive Groups from RDChem Mol.")

        try:
           molecule.eleCompFromMol()
        except:
            window.showErrorMessage("Determining Elemental Compostion from RDChem Mol.")
        
        try:
           molecule.CHOFromEleComp()
        except:
            window.showErrorMessage("Determining CHO from Elemental Compostion.")

        try:        
           molecule.RoSFromEleComp()
        except:
            window.showErrorMessage("Calculating Rule of Six from Elemental Compostion.")

        try:
           molecule.OBFromEleComp()
        except:
            window.showErrorMessage("Calculating Oxygen Balance from Elemental Compostion.")

        try:
           molecule.oreoOnsetTadjustment()
        except:
            window.showErrorMessage("Adjusting O.R.E.O.S. Calculation for Onset Temperature")

        try:
           molecule.oreoSafeScaleCal()
        except:
            window.showErrorMessage("Determining O.R.E.O.S. Hazard by Scale.")

        try:
           molecule.Td24FromThermalProps()
        except:
            window.showErrorMessage("Calculating T<sub>D24</sub> from Thermal Properties.")

        try:
           molecule.yoshidaFromThermalProps()
        except:
            window.showErrorMessage("Calculating Yoshida values from Thermal Properties.")

        try:
            molecule.makeFolderForMolData()
        except:
            window.showErrorMessage("Making Folder to Hold Molecule Data Files.")

    def createReport(self):
        #try:
        moleculeData = self.render_molecule()
           #img = moleculeData.molToIMG()
           #results = asdict(moleculeData)
           #create_pdf(results["name"], results, img) #results["molIMG"]) 
        imageData = moleculeData.molToBytes()
        dataURL = 'data:image/png;base64,' + imageData
        mdReportCreation(moleculeData, dataURL)
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)
        #except:
            #window.showErrorMessage("Generating Memo PDF from given values.")

    def writeToDatabase(self, molecule, Database):
        #molecule.genAdditionalValues()
        selectedMolData = cleanMolDataFrame(molecule)
        storedData = pd.read_csv(Database, index_col=0)
        checkData = pd.read_csv(Database)
        print('\n\n\n')
        if selectedMolData['SMILES'][0] in storedData.index:
            print('found')
            userInteract = self.interactiveErrorMessage('Molecule Already in Database. Would you like to overwrite it?')
            if userInteract == QMessageBox.Yes:
                storedData.update(selectedMolData)
                outputData = storedData
                row_index = checkData.index[checkData['SMILES'] == selectedMolData['SMILES'][0]].tolist()
                print(row_index)
                folderInfo = checkData['dataFolder'][row_index[0]]
                print(folderInfo)
                output_index = outputData.index[checkData['SMILES'] == selectedMolData['SMILES'][0]].tolist()
                outputData['dataFolder'][output_index[0]] = folderInfo

            elif userInteract == QMessageBox.No:
                outputData = storedData

            else:
                outputData = storedData

        else:
            outputData = pd.concat([storedData, selectedMolData])

        outputData['SMILES'] = outputData.index
        outputData = outputData[ ['SMILES'] + [ col for col in outputData.columns if col != 'SMILES' ] ]
        print(outputData)
        outputData.to_csv(Database, index=False)

    def getColorForValue(self, hazardClass):
        # Color-coding logic
        if hazardClass == 'High Hazard':
            return QColor(255, 0, 0)  # Red
        elif hazardClass == 'Medium Hazard':
            return QColor(255, 255, 0)  # Yellow
        elif hazardClass == 'Low Hazard':
            return QColor(0, 255, 0)  # Green
        else:
            return QColor(0, 0, 255)  # Blue

    def genMoleculeFromUserInput(self):
        # Get the SMILES string from the input field
        smiles = self.smiles_input.text()
        name = self.name_input.text()
        mp = self.mp_input.text()
        mpEnd = self.mpEnd_input.text()
        Qdsc = self.Qdsc_input.text()
        QUnits = self.QUnitsSelection.currentText()
        TE = self.TE_input.text()
        Tinit = self.Tinit_input.text()
        proj = self.proj_input.text()

        dataFolder = self.findFolder(defaultDB)

        # Create an RDKit molecule from the SMILES string
        addedMolecule = thermalDexMolecule(SMILES=smiles, name=name, mp=mp, mpEnd=mpEnd, Q_dsc=Qdsc, Qunits=QUnits, onsetT=TE, initT=Tinit, proj=proj, dataFolder=dataFolder)
        return addedMolecule

    def checkIfSMILESAreValid(self, molecule):
        if molecule.SMILES != '' and molecule.SMILES is not None:
            molecule.mol = MolFromSmiles(molecule.SMILES)

        else:
            molecule.mol = None

    def displayTheMolecule(self, molecule, display):
        # Make Pixmap Image to Display.
        pixmap = molecule.molToQPixmap()
        scene = QGraphicsScene()
        scene.addPixmap(pixmap)
        display.setScene(scene)
        molecule.molPixmap = None

    def render_molecule(self):
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)

        # Generate a thermalDexMolecule from the user input.
        addedMolecule = self.genMoleculeFromUserInput()

        # Check that the user provided SMILES are valid.
        self.checkIfSMILESAreValid(addedMolecule)
        
        if addedMolecule.mol is not None:
            # Make Pixmap Image to Display.
            self.displayTheMolecule(addedMolecule, self.mol_display)

            # Calculate Core Properties
            self.genCoreValuesFromMol(addedMolecule)

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
                self.Td24Label.setText('T<sub>D24</sub>: <b>' + d24Str + ' °C</b>')
            if addedMolecule.onsetT != 'nan' and addedMolecule.onsetT != '' and addedMolecule.onsetT != None and addedMolecule.onsetT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Onset Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 201
            if addedMolecule.initT != 'nan' and addedMolecule.initT != '' and addedMolecule.initT != None and addedMolecule.initT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Initiation Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 202

            #createDatabase(addedMolecule)
            self.writeToDatabase(addedMolecule, defaultDB)
            print('\n')
            print(addedMolecule)
            print('\n')
            fileCounter = self.countFiles(defaultDB)
            self.filesCount.setText(f"{str(fileCounter)} Attached Files")
            self.attachedFilesLabel.show()
            self.filesCount.show()
            self.attach_button.show()
            tempFiles = [f for f in listdir('./_core/UserAddedData/temp') if isfile(join('./_core/UserAddedData/temp', f))]
            for file in tempFiles:
                try:
                    remove(f'./_core/UserAddedData/temp/{file}')
                except:
                    pass
            return addedMolecule


        else:
            # Display an error message if the SMILES string is invalid
            self.error_message = QLabel('Invalid SMILES string. Please enter a valid SMILES.')
            #layout = self.layout()
            layout.addWidget(self.error_message)
            self.error_flag = 100

if __name__ == '__main__':
    defaultDB, highEnergyGroups, expEnergyGroups, yosidaMethod = altImportConfig()
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('.\\_core\\ThermalDexIcon.ico'))
    window = MolDrawer()
    window.show()
    window.raise_()
    window.activateWindow()
    sys.exit(app.exec_())
