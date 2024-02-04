
from markdown_pdf import MarkdownPdf, Section
#from custom_Markdown_PDF import MarkdownPdf, Section
from thermDexMolecule import thermalDexMolecule

highEnergyGroups = ['C', 'N=N']
expEnergyGroups = ['N', 'OO']

def mdReportCreation(molecule, dataURL):
    memo = MarkdownPdf(toc_level=1)

    #<style>
    css = '''
    img {
    display: block;
    margin-left: auto;
    margin-right: auto;
    }
    '''
    #</style>

    memo.add_section(Section("# Thermal Hazard Assessment Memo\n", toc=False)) #, css)

    report_body = f'''
## {molecule.name}
![../ThermalDexSplash.jpg](../ThermalDexSplash.jpg)
                 
## Molecule Properties 
SMILES: {molecule.SMILES}\n
Q<sub>DSC</sub> = 550 J g-1
T<sub>D24</sub> = 130 °C

<div style="text-align: center;">
<img src="{dataURL}" alt="HTML image test" width="50%"/>
</div>

## Interpritation 
These results have been calculated using X and they show Y.
'''

    memo.add_section(Section(report_body)) #, css)
    memo.save("altmemo.pdf")

molecule = thermalDexMolecule(SMILES='c1ccc(C)cc1',name='Test mol Toluene')
molecule.genMol()
molecule.molToIMG() #.genAllValues() #highEnergyGroups,expEnergyGroups)
imageData = molecule.molToBytes()
dataURL = 'data:image/png;base64,' + imageData
print(dataURL)
mdReportCreation(molecule, dataURL)
