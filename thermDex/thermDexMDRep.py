
from markdown_pdf import MarkdownPdf, Section
#from custom_Markdown_PDF import MarkdownPdf, Section
import fitz
from thermDexMolecule import thermalDexMolecule

highEnergyGroups = ['C', 'N=N']
expEnergyGroups = ['N', 'OO']

def mdReportCreation(molecule, dataURL):
    memo = MarkdownPdf(toc_level=0)

    #<style>
    css = '''
    img {
    display: block;
    margin-left: auto;
    margin-right: auto;
    }
    '''
    #</style>

    #memo.add_section(Section("# Thermal Hazard Assessment Memo\n", toc=False)) #, css)

    if molecule.mp == '' and molecule.mpEnd == '':
        mpString = ''

    elif molecule.mp != '' and molecule.mpEnd == '':
        mpString = f'mp: {molecule.mp} °C'

    elif molecule.mp != '' and molecule.mpEnd != '':
        mpString = f'mp: {molecule.mp} to {molecule.mpEnd} °C'

    report_body = f'''
<div style="text-align: center;">

# Thermal Hazard Assessment Memo\n

</div>

## {molecule.name}
<div style="text-align: center;">
<img src="{dataURL}" alt="HTML image test" width="50%"/>
</div>
                 
### Molecule Properties 
SMILES: {molecule.SMILES}\n
Formula: {molecule.eleComp}\n
{mpString}\n

### Results
<table align="Center" style="border-spacing: 10px;">
    <tr>
        <td colspan="3">High Energy Groups: ({molecule.HEG}) {molecule.HEG_list}</td>
    </tr>
    <tr>
        <td colspan="3">Explosive Groups: ({molecule.EFG}) {molecule.EFG_list}</td>
    </tr>
    <tr>
        <td>Rule of Six = {molecule.RoS_val}</td>
        <td>Oxygen Balance = {molecule.OB_val}</td>
        <td> </td>
    </tr>
    <tr>
        <td>Q<sub>DSC</sub> = {molecule.Q_dsc} {molecule.Qunits}</td>
        <td>T<sub>onset</sub> = {molecule.onsetT}</td>
        <td>T<sub>init</sub> = {molecule.initT}</td>
    </tr>
    <tr>
        <td>Impact Sensitivity = {molecule.IS_val}</td>
        <td>Explosive Propagation = {molecule.EP_val}</td>
        <td>T<sub>D24</sub> = {molecule.Td24} °C</td>
    </tr>
</table>

<br>

<table align="Center" style="border-collapse:collapse; border:1px solid black; "border-spacing:20px;">
    <tr style="border:1px solid black; border-collapse: collapse;">
        <td style="border:1px solid black; border-collapse: collapse;"> <5 g</td>
        <td style="border:1px solid black; border-collapse: collapse;"> 5 to 100 g</td>
        <td style="border:1px solid black; border-collapse: collapse;"> 100 to 500 g</td>
        <td style="border:1px solid black; border-collapse: collapse;"> >500 g</td>
    </tr>
    <tr style="border:1px solid black; border-collapse: collapse; padding:0px;">
        <td style="border:1px solid black; border-collapse: collapse;">{molecule.oreoSmallScale_des}</td>
        <td style="border:1px solid black; border-collapse: collapse;">{molecule.oreoTensScale_des}</td>
        <td style="border:1px solid black; border-collapse: collapse;">{molecule.oreoHundsScale_des}</td>
        <td style="border:1px solid black; border-collapse: collapse;">{molecule.oreoLargeScale_des}</td>
    </tr>
</table>

### Interpretation
These results have been calculated using X<sup>1</sup> and they show Y<sup>2</sup>.

<small>[1]: *Org. Proc. Res. Dev.,* 2011, 2341-2356</small>\n
<small>[2]: *Org. Proc. Res. Dev.,* 2011, 2117-2119</small>
'''

    memo.add_section(Section(report_body, toc=False)) #, css)
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
