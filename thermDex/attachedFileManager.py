import sys
import shutil
import subprocess
from os import listdir, remove, system
from os.path import isfile, join
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QListWidget, QPushButton, QLabel, QFileDialog, QMessageBox

class FileManagementWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Attached Files")
        self.setGeometry(100, 100, 400, 300)
        self.fileDir = './_core/UserAddedData'

        self.file_list_widget = QListWidget()
        self.load_files()

        self.add_button = QPushButton("+ Add File")
        self.remove_button = QPushButton("- Remove File(s)")
        self.open_button = QPushButton("Open File(s)")
        
        self.add_button.clicked.connect(self.add_file)
        self.remove_button.clicked.connect(self.remove_file)
        self.open_button.clicked.connect(self.open_file)

        layout = QVBoxLayout()
        layout.addWidget(QLabel("Files in the Folder:"))
        layout.addWidget(self.file_list_widget)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.add_button)
        button_layout.addWidget(self.remove_button)
        open_button_layout = QHBoxLayout()
        open_button_layout.addWidget(self.open_button)
        layout.addLayout(button_layout)
        layout.addLayout(open_button_layout)
        self.setLayout(layout)

    def open_window(self):
        self.newWindow = FileManagementWindow(self)
        self.newWindow.show()

    def load_files(self):
        # Simulating loading files from a folder
        #files = ["file1.txt", "file2.txt", "file3.txt"]
        files = [f for f in listdir(self.fileDir) if isfile(join(self.fileDir, f))]
        self.file_list_widget.addItems(files)

    def add_file(self):
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFiles)
        if file_dialog.exec_():
            selected_files = file_dialog.selectedFiles()
            for file in selected_files:
                shutil.copy2(file, self.fileDir)
            self.file_list_widget.clear()
            self.load_files()

    def remove_file(self):
        deleteMsg = QMessageBox()
        deleteMsg.setIcon(QMessageBox.Information)
        deleteMsg.setText("Action Needed")
        deleteMsg.setInformativeText('Are you sure you want to delete selected files?')
        deleteMsg.setWindowTitle("ThermalDex - Delete Files")
        deleteMsg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        userConfirm = deleteMsg.exec()
        if userConfirm == QMessageBox.Yes:
            for item in self.file_list_widget.selectedItems():
                self.file_list_widget.takeItem(self.file_list_widget.row(item))
                itemPath = self.fileDir + '/' + item.text()
                remove(itemPath)

        elif userConfirm == QMessageBox.No:
            pass

        else:
            pass


    def open_file(self):
        for item in self.file_list_widget.selectedItems():
            print(item.text())
            itemPath = self.fileDir + '/' + item.text()
            subprocess.run(['start', itemPath], check=True, shell=True)

