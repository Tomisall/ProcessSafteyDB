from markdown_pdf import MarkdownPdf, Section
#from custom_Markdown_PDF import MarkdownPdf, Section

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

report_body = '''## Name
![../ThermalDexSplash.jpg](../ThermalDexSplash.jpg)
                 
## Molecule Properties 
Q<sub>DSC</sub> = 550 J g-1
T<sub>D24</sub> = 130 °C

<div style="text-align: center;">
<img src="../ThermalDexSplash.jpg" alt="HTML image test" width="50%"/>
</div>

## Interpritation 
These results have been calculated using X and they show Y.
'''

memo.add_section(Section(report_body)) #, css)
memo.save("altmemo.pdf")
