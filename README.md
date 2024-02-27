# ThermalDex
[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
![Qt](https://img.shields.io/badge/Qt-%23217346.svg?style=for-the-badge&logo=Qt&logoColor=white)
![HotCookie](./HouseKeeping/GitHubBadges/HotCookie.svg)

ThermalDex is an application developed in Python using RDKit and PyQt5 that facilitates the identification and assessment of thermal hazards in molecules. It offers suite of tools for analyzing chemical structures, calculating safety metrics, and generating detailed reports to aid in the safe handling and processing of potentially hazardous materials.

## Features

- **SMILES Input**: ThermalDex accepts input in the form of SMILES (Simplified Molecular Input Line Entry System), allowing users to input molecular structures easily. ThermalDex will also accept direct copy/paste from Chemaxon's Marvin/JChem files (ChemDraw Support is *ToDo*).

- **High Energy and Explosive Group Identification**: The application identifies high-energy and potentially explosive groups within molecules, providing valuable insights into their thermal stability and safety. This is achived via a substructure search and a .csv of the relavent groups.

- **O.R.E.O.S. Assessment**: ThermalDex calculates O.R.E.O.S. assessment values, incorporating safety metrics such as Oxygen Balance and Rule of Six to evaluate the potential hazards associated with a molecule.

- **DSC Data Analysis**: The application can record and analyze the results from Differential Scanning Calorimetry (DSC), enabling users to determine critical parameters such as T<sub>D24</sub> (Temperature of Decomposition), Yoshida Impact Sensitivity, and Yoshida Explosive Propagation.

- **Graphical Visualization**: ThermalDex features a Qt GUI that displays molecular structures graphically, allowing users to visualize the chemical composition and hazard potential of each molecule. Hazard levels color-coded for intuitive interpretation.

- **PDF Report Generation**: The application generates detailed PDF reports of thermal hazard assessments for each molecule, providing comprehensive documentation of the analysis results.

- **Local Storage**: ThermalDex maintains a user-friendly local .csv file 'database' (cough) that stores all assessment results, making it easy to track and manage data over time.

## Getting Started

### Executable Releases 
Exicutables will be Released via Github Releases. (Currently we are testing internally but will release in near future). These executables have be generated thanks to Pyinstaller. ThermalDex has been a Windows first project due to the nature of the devices used for testing, but I will get a MacOS exicutable shortly and aim to generate at least a .deb package for Linux down the line. If compatabliity issues are encountered I would recomemend looking to running in Python. 

### Python Code

To get started with ThermalDex, follow these steps:

1. Clone the repository to your local machine:
   ```
   git clone https://github.com/Tomisall/ThermalDex.git
   ```

2. Install the required dependencies using pip:
   ```
   pip install -r requirements.txt
   ```

3. Launch the application:
   ```
   python main.py
   ```

4. Input the SMILES notation of the molecule you wish to analyze and explore the various features and tools offered by ThermalDex.

## Contributing

Contributions to ThermalDex are welcome! If you encounter any issues, have suggestions for improvements, or would like to contribute new features, please feel free to submit a pull request or open an issue on GitHub.

Please get in touch if you have any issues or would like to see any additional features. I want ThermalDex to be a really useful tool for those thinking about (and those who don't want to but should think about) Thermal Hazards in organic synthesis.

## License

ThermalDex is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.