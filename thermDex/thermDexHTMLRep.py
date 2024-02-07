import sys
import subprocess
from pathlib import Path
from xhtml2pdf import pisa
import time

def convert_html_to_pdf(html_string, pdf_path):
    with open(pdf_path, "wb") as pdf_file:
        pisa_status = pisa.CreatePDF(html_string, dest=pdf_file)
        
    return not pisa_status.err

def mdReportCreation(molecule, dataURL):
    if molecule.mp == '' and molecule.mpEnd == '':
        mpString = ''

    elif molecule.mp != '' and molecule.mpEnd == '':
        mpString = f'mp: {molecule.mp} °C'

    elif molecule.mp != '' and molecule.mpEnd != '':
        mpString = f'mp: {molecule.mp} to {molecule.mpEnd} °C'

    html_content = f'''
<!DOCTYPE html>
<html>
<head>
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=UTF-8">
</head>
<body>
<p><link rel="stylesheet" href="./_core/style.css"></p>
<div style="text-align: center;">
<h1>Thermal Hazard Assessment Memo</h1>
</div>
<h2>{molecule.name}</h2>
<div style="text-align: center;">
<img src="{dataURL}" alt="HTML image test" width="50%"/>
</div>              
<h3>Molecule Properties</h3>
<p>SMILES: {molecule.SMILES}</p>
<p>Formula: {molecule.eleComp}</p>
<p>{mpString}</p>
<h3>Results</h3>

<table align="Center" style="border-spacing: 10px; width: 100%;">
<tr class="even">
<td colspan="3">High Energy Groups: ({molecule.HEG}) {molecule.HEG_list} &nbsp;</td>
</tr>
<tr>
<td colspan="3">Explosive Groups: ({molecule.EFG}) {molecule.EFG_list}</td>
</tr>
<tr class="even">
<td>Rule of Six = {molecule.RoS_val}</td>
<td>Oxygen Balance = {molecule.OB_val}</td>
<td> </td>
</tr>
<tr class="even">
<td>Q<sub>DSC</sub> = {molecule.Q_dsc} {molecule.Qunits}</td>
<td>T<sub>onset</sub> = {molecule.onsetT}</td>
<td>T<sub>init</sub> = {molecule.initT}</td>
</tr>
<tr class="even">
<td>Impact Sensitivity = {molecule.IS_val}</td>
<td>Explosive Propagation = {molecule.EP_val}</td>
<td>T<sub>D24</sub> = {molecule.Td24} °C</td>
</tr>
</table>

<p>O.R.E.O.S. assessment of risk by scale:</p>
<table>
<thead>
<tr class="header">
<th> <5 g</th>
<th> 5 to 100 g</th>
<th> 100 to 500 g</th>
<th> >500 g</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>{molecule.oreoSmallScale_des}</td>
<td>{molecule.oreoTensScale_des}</td>
<td>{molecule.oreoHundsScale_des}</td>
<td>{molecule.oreoLargeScale_des}</td>
</tr>
</tbody>
</table>
<h3>Interpretation</h3>
<p>These results have been calculated using X<sup>1</sup> and they show Y<sup>2</sup>.</p>
<p><small>[1]: <i>Org. Proc. Res. Dev.,</i> 2011, 2341-2356</small><br>
<small>[2]: <i>Org. Proc. Res. Dev.,</i> 2011, 2117-2119</small></p>
</body>
'''
    convert_html_to_pdf(html_content, "./AssessmentMemos/htmlMemo.pdf")
    time.sleep(2)
    subprocess.run(['start', "./AssessmentMemos/htmlMemo.pdf"], check=True, shell=True)