from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph

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


def create_pdf(Name, results, interpretation, image_path):
    filename = "thermal_assessment_report.pdf"

    # Create the PDF
    c = canvas.Canvas(filename, pagesize=A4)

    # Add title
    c.setFont("Helvetica", 16)
    c.drawCentredString(300, 775, "Thermal Hazard Assessment Report")
    

    # Name Molecule
    c.setFont("Helvetica", 14)
    c.drawString(50, 725, Name)

    # Add image at the top center
    if image_path:
        c.drawInlineImage(image_path, 150, 405, width=300, height=300)

    # Add results section
    c.setFont("Helvetica", 12)
    c.drawString(50, 370, "Results:")

    # Add the actual results from your assessment
    for i, (key, value) in enumerate(results.items()):
        c.drawString(70, 350 - i * 20, f"{key}: {value}")

    # Add interpretation section
    c.drawString(50, 200, "Interpretation:")
    c.setFont("Helvetica", 12)
    # Add your interpretation text
    interpretation_text = interpretation
    # interpPara = Paragraph(interpretation_text, my_Style)
    # interpPara.drawOn(c, 70, 480)

    i = 180
    for line in interpretation_text:
        c.drawString(70, i, line)
        i -= 20

    # Add footer with disclaimer
    c.setFont("Helvetica", 8)
    disclaimer_text = "This report may contain confidential information."
    c.drawCentredString(300, 30, disclaimer_text)

    # Save the PDF
    c.save()
    print(f"PDF report generated: {filename}")

# Example usage:
Name = "Dave"
results = {"Parameter1": "Value1", "Parameter2": "Value2", "Parameter3": "Value3"}
interpretation = ["This is my assessment of the molecule", "On balance:", "Confrimed to be deadly"]
image_path = "./_core/ThermalDexIcon.jpg"

create_pdf(Name, results, interpretation, image_path)
