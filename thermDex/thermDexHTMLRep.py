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
<div class="imgDiv">
<img src="{dataURL}" alt="HTML image test" style="object-fit: cover;"/>
</div>              
<h3>Molecule Properties</h3>
<p>SMILES: {molecule.SMILES}<br>
Formula: {molecule.eleComp}<br>
MW: {"{:.2f}".format(molecule.MW)} g mol<sup>-1</sup><br>
{mpString}</p>
<h3>Results</h3>
<table align="Center" class="secretTable">
<tr class="secretTable">
<td colspan="3" class="secretTable">High Energy Groups: ({molecule.HEG}) {', '.join(molecule.HEG_list)} &nbsp;</td>
</tr>
<tr class="secretTable">
<td colspan="3" class="secretTable">Explosive Groups: ({molecule.EFG}) {', '.join(molecule.EFG_list)}</td>
</tr>
<tr class="secretTable" class="secretTable">
<td class="secretTable">Rule of Six = {molecule.RoS_val}</td>
<td class="secretTableMid">Oxygen Balance = {"{:.2f}".format(molecule.OB_val)}</td>
<td class="secretTable"> </td>
</tr>
<tr class="secretTable">
<td class="secretTable">Q<sub>DSC</sub> = {molecule.Q_dsc} {molecule.Qunits[:-2]}<sup>-1</sup></td>
<td class="secretTableMid">T<sub>onset</sub> = {molecule.onsetT} °C</td>
<td class="secretTable">T<sub>init</sub> = {molecule.initT} °C</td>
</tr>
<tr class="secretTable">
<td class="secretTable">Impact Sensitivity = {"{:.2f}".format(molecule.IS_val)}</td>
<td class="secretTableMid">Explosive Propagation = {"{:.2f}".format(molecule.EP_val)}</td>
<td class="secretTable">T<sub>D24</sub> = {molecule.Td24} °C</td>
</tr>
</table>

<p>O.R.E.O.S. assessment of risk by scale:</p>
<table align="Center" width="65%" >
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
<div id="footer_content" style="text-align: center;">
<p><small>This memo may contain confidential information.</small><br>
<small><a href="https://github.com/Tomisall/ProcessSafteyDB">ThermalDex</a></small></p>
</div>
</body>
'''
    convert_html_to_pdf(html_content, "./AssessmentMemos/htmlMemo.pdf")
    time.sleep(2)
    subprocess.run(['start', "./AssessmentMemos/htmlMemo.pdf"], check=True, shell=True)