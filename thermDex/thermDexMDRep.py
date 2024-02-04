
from markdown_pdf import MarkdownPdf, Section
#from custom_Markdown_PDF import MarkdownPdf, Section
import fitz
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

    doc = fitz.open("altmemo.pdf")
    numpages = doc.page_count  # number of pages
    footer_text = "Page %i of %i"
    header_text = "My Matplotlib File"
    blue = fitz.pdfcolor["blue"]

    for page in doc:
        prect = page.rect
        header_rect = fitz.Rect(0, 36, prect.width, 56)  # height 20 points
        page.insert_textbox(header_rect, header_text,
                            fontname="hebo", color=blue,
                            align=fitz.TEXT_ALIGN_CENTER)

        ftext = footer_text % (page.number + 1, numpages)
        y1 = prect.height - 36  # bottom of footer rect
        y0 = y1 - 20  # top of footer rect
        footer_rect = fitz.Rect(0, y0, prect.width, y1)  # rect has full page width
        page.insert_textbox(footer_rect, 'This report may contain confidential information', align=fitz.TEXT_ALIGN_CENTER)

        doc.save("altmemo_2.pdf")

molecule = thermalDexMolecule(SMILES='c1ccc(C)cc1',name='Test mol Toluene')
molecule.genMol()
molecule.molToIMG() #.genAllValues() #highEnergyGroups,expEnergyGroups)
imageData = molecule.molToBytes()
dataURL = 'data:image/png;base64,' + imageData
print(dataURL)
mdReportCreation(molecule, dataURL)
