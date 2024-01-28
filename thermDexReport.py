from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from thermDex.thermDexMolecule import thermalDexMolecule
from dataclasses import asdict

#arial_path = "/System/Library/Fonts/Supplemental/Arial.ttf"
#arial = TTFont("Arial","/System/Library/Fonts/Supplemental/Arial.ttf")
pdfmetrics.registerFont(TTFont('Arial', './_core/Arial.ttf'))
memoFont = "Arial"

my_Style=ParagraphStyle('My Para style',
fontName='Times-Roman',
backColor='#F1F1F1',
fontSize=16,
borderColor='#FFFF00',
borderWidth=2,
borderPadding=(20,20,20),
leading=20,
alignment=0
)

def apply_scripting(textobject, text, rise):
    textobject.setFont(memoFont, 8)
    textobject.setRise(rise)
    textobject.textOut(text)
    textobject.setFont(memoFont, 11)
    textobject.setRise(0)


def create_pdf(Name, results, interpretation, image_path):
    filename = "thermal_assessment_report.pdf"

    # Create the PDF
    c = canvas.Canvas(filename, pagesize=A4)

    fontsAva = c.getAvailableFonts()
    print(fontsAva)

    # Add title
    c.setFont(memoFont, 16)
    c.drawCentredString(300, 775, "Thermal Hazard Assessment Report")
    
    # Name Molecule
    c.setFont(memoFont, 14)
    c.drawString(50, 725, Name)

    # Add image at the top center
    if image_path:
        c.drawInlineImage(image_path, 150, 490, width=225, height=225)

    # Add properties section
    c.setFont(memoFont, 11)
    c.drawString(50, 455, "Properties:")

    # Add the actual properties from your assessment
    c.drawString(70, 435, "SMILES: " + results["SMILES"])
    try:
        c.drawString(70, 415, "Name: " + results["name"])
    except:
        c.drawString(70, 415, "Name: ")
    try:
        c.drawString(70, 395, "Formular: " + results["molForm"])
    except:
        c.drawString(70, 395, "Formular: ")   
    try:
        c.drawString(70, 375, "mp: " + str(results["mp"]) + " to " +  str(results["mpEnd"]))
        textobject = c.beginText()
        textobject.setTextOrigin(200, 375)
        textobject.textOut('37')
        apply_scripting(textobject, '%', -4)
        c.drawText(textobject)
    except:
        c.drawString(70, 375, "mp: " + " " + " to " + " ")

    # Add results section
    c.setFont(memoFont, 11)
    c.drawString(50, 335, "Results:")

    # Add the actual results from your assessment
    # for i, (key, value) in enumerate(results.items()):
    #    c.drawString(70, 315 - i * 20, f"{key}: {value}")
    c.drawString(70, 315, "High Energy Groups =  " + str(results["HEG"]) + " (" + ", ".join(results["HEG_list"]) + ")")

    textobject = c.beginText()
    textobject.setTextOrigin(70, 295)
    textobject.textOut("Q")   
    apply_scripting(textobject, "DSC", -4)
    textobject.textOut(" = " + str(results["Q_dsc"]))
    c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(250, 295)
    textobject.textOut("T")   
    apply_scripting(textobject, "onset", -4)
    textobject.textOut(" = " + str(results["onsetT"]))
    c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(430, 295)
    textobject.textOut("T")   
    apply_scripting(textobject, "init", -4)
    textobject.textOut(" = " + str(results["initT"]))
    c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(70, 275)
    textobject.textOut("Rule of Six")   
    #apply_scripting(textobject, "init", -4)
    textobject.textOut(" = " + str(results["RoS_val"]))
    c.drawText(textobject)

    textobject = c.beginText()
    textobject.setTextOrigin(250, 275)
    textobject.textOut("Oxygen Balance")   
    #apply_scripting(textobject, "init", -4)
    textobject.textOut(" = " + str(results["OB_val"]))
    c.drawText(textobject)



    # Add interpretation section
    c.drawString(50, 200, "Interpretation:")
    c.setFont(memoFont, 12)
    # Add your interpretation text
    interpretation_text = interpretation
    # interpPara = Paragraph(interpretation_text, my_Style)
    # interpPara.drawOn(c, 70, 480)

    i = 180
    for line in interpretation_text:
        c.drawString(70, i, line)
        i -= 20

    # Add footer with disclaimer
    c.setFont(memoFont, 8)
    disclaimer_text = "This report may contain confidential information."
    c.drawCentredString(300, 30, disclaimer_text)

    # Save the PDF
    c.save()
    print(f"PDF report generated: {filename}")

# Example usage:
TNT = thermalDexMolecule(SMILES='CC1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', name='TNT', Q_dsc=770.26, onsetT=90.14, initT=74.62, proj='PDF_Test')
TNT.genAllValues()
Name = TNT.name
results = asdict(TNT) #{"Parameter1": "Value1", "Parameter2": "Value2", "Parameter3": "Value3"}
interpretation = ["This is my assessment of the molecule", "On balance:", "Confrimed to be deadly"]
image_path = TNT.molPixmap #"./_core/ThermalDexIcon.jpg"

create_pdf(Name, results, interpretation, image_path)
