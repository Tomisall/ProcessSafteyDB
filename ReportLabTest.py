from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

def create_pdf(results, interpretation):
    filename = "thermal_assessment_report.pdf"

    # Create the PDF
    c = canvas.Canvas(filename, pagesize=letter)

    # Add title
    c.setFont("Helvetica", 16)
    c.drawCentredString(300, 750, "Thermal Hazard Assessment Report")

    # Add results section
    c.setFont("Helvetica", 12)
    c.drawString(50, 720, "Results:")

    # Add the actual results from your assessment
    for i, (key, value) in enumerate(results.items()):
        c.drawString(70, 700 - i * 20, f"{key}: {value}")

    # Add interpretation section
    c.drawString(50, 600, "Interpretation:")
    c.setFont("Helvetica", 10)
    # Add your interpretation text
    interpretation_text = (
        "This section provides a simple interpretation of the results."
        " You can customize this based on your specific assessment."
    )
    c.drawParagraph(interpretation_text, 70, 580)

    # Add footer with disclaimer
    c.setFont("Helvetica", 8)
    disclaimer_text = "This report may contain confidential information."
    c.drawCentredString(300, 30, disclaimer_text)

    # Save the PDF
    c.save()
    print(f"PDF report generated: {filename}")

# Example usage:
results = {"Parameter1": "Value1", "Parameter2": "Value2", "Parameter3": "Value3"}
interpretation = "This is a sample interpretation text."

create_pdf(results, interpretation)