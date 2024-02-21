import sys
import subprocess
from pathlib import Path
from xhtml2pdf import pisa
import time
from os import listdir, remove
from os.path import isfile, join, abspath
import fitz
from PIL import Image
from io import BytesIO
import base64
from docx2pdf import convert
from datetime import datetime
import win32com.client
import pythoncom

def convert_html_to_pdf(html_string, pdf_path):
    with open(pdf_path, "wb") as pdf_file:
        pisa_status = pisa.CreatePDF(html_string, dest=pdf_file)
        
    return not pisa_status.err

def mdReportCreation(molecule, dataURL):
    wdFormatPDF = 17

    if molecule.mp == '' and molecule.mpEnd == '':
        mpString = ''

    elif molecule.mp != '' and molecule.mpEnd == '':
        mpString = f'mp: {molecule.mp} °C'

    elif molecule.mp != '' and molecule.mpEnd != '':
        mpString = f'mp: {molecule.mp} to {molecule.mpEnd} °C'


    additionalData = '<h2>Additional Data</h2>'
    attachedFiles = [f for f in listdir(molecule.dataFolder) if isfile(join(molecule.dataFolder, f))]
    print('\n\n')
    print(molecule.dataFolder)
    print(attachedFiles)
    imgList = ['.png', '.PNG', '.jpeg', '.JPEG', '.jpg', '.JPG', '.svg', '.SVG']
    pdfList = ['.pdf', '.PDF']
    docList = ['.doc', '.docx']
    for file in attachedFiles:
        if file.endswith(tuple(imgList)):
             print(file)
             additionalData += f'''
<h3>{file}</h3>
<div class="attachedFile">
<img src="{molecule.dataFolder}/{file}" alt="{file}" style="object-fit: cover;"/>
</div>
'''
        elif file.endswith(tuple(pdfList)):
            additionalData += f'<h3>{file}</h3>'
            pdfFile = fitz.open(f'{molecule.dataFolder}/{file}')
            noOfPages = pdfFile.page_count
            for pageNo in range(noOfPages):
                page = pdfFile.load_page(pageNo) 
                pix = page.get_pixmap()
                #image_data = pix.tobytes()
                #byte_array = BytesIO()
                #byte_array.write(image_data, format='PNG')
                image_bytes = pix.tobytes("PNG")
                pngOfPage = base64.b64encode(image_bytes).decode('utf-8')
                pageDataURL = 'data:image/png;base64,' + pngOfPage
                additionalData += f'''
<div class="attachedFile">
<img src="{pageDataURL}" alt="{file}" style="object-fit: cover;"/>
</div>
'''

        elif file.endswith(tuple(docList)):
            additionalData += f'<h3>{file}</h3>'
            now = datetime.now()
            neatNow = now.strftime("%d-%b-%Y_%H-%M-%S")
            outName = f'./_core/UserAddedData/temp/temp_{str(neatNow)}.pdf'
            #convert(f'{molecule.dataFolder}/{file}', outName)
            absolutePath = abspath(f'{molecule.dataFolder}/{file}')
            outAbs = abspath(outName)
            word = win32com.client.Dispatch('Word.Application', pythoncom.CoInitialize())
            doc = word.Documents.Open(absolutePath)
            doc.SaveAs(outAbs, FileFormat=wdFormatPDF)
            doc.Close()
            #word.Quit()
            pdfFile = fitz.open(outName)
            noOfPages = pdfFile.page_count
            for pageNo in range(noOfPages):
                page = pdfFile.load_page(pageNo) 
                pix = page.get_pixmap()
                #image_data = pix.tobytes()
                #byte_array = BytesIO()
                #byte_array.write(image_data, format='PNG')
                image_bytes = pix.tobytes("PNG")
                pngOfPage = base64.b64encode(image_bytes).decode('utf-8')
                pageDataURL = 'data:image/png;base64,' + pngOfPage
                additionalData += f'''
<div class="attachedFile">
<img src="{pageDataURL}" alt="{file}" style="object-fit: cover;"/>
</div>
'''
        else:
            print(f'poo poo {file}')



    Td24Formated = "{:.2f}".format(molecule.Td24) + " °C"
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
<td class="secretTable">T<sub>D24</sub> = {Td24Formated}</td>
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
    if additionalData != '<h2>Additional Data</h2>':
        html_content += additionalData
    else:
        print('\n\nFoo\n\n')
    convert_html_to_pdf(html_content, "./AssessmentMemos/htmlMemo.pdf")
    time.sleep(1)
    subprocess.run(['start', "./AssessmentMemos/htmlMemo.pdf"], check=True, shell=True)