from fpdf import FPDF
import MolDisplay
import os
import tempfile
import numpy as np
import pandas as pd

class PDF_generator (FPDF):
    filepath = './test.pdf'
    hdrimg = ''
    hdrtxt = ''
    def __init__(self, filepath = './_test.pdf', hdrtxt=None, hdrimg=None):
        super().__init__()
        self.hdrtxt = hdrtxt
        self.hdrimg = hdrimg
        self.filepath=filepath

    def init_pdf (self):
        self.set_author('CIxTools')
        self.add_page()
    def save (self):
        self.output (self.filepath)
    def add_text (self, txtval, align='L', color=(0,0,0), bold='', fontsize = 16):
        self.set_font('Arial', bold, fontsize) #B is bold
        self.set_text_color(color[0], color[1], color[2])
        self.multi_cell(0, 10, txt=txtval, align=align)
    def add_image_fromfile (self, imgpath):
        self.image(imgpath, link='', type='')
    def add_MoleculeFromSMILES (self, smileslist):
        temp_file_path = tempfile.NamedTemporaryFile(suffix=".png").name
        MolDisplay.ShowMols(smileslist, fname=temp_file_path, subImgSize=(400,400))
        self.image(temp_file_path, w=200)
        os.remove (temp_file_path)
    def border(self):
        self.set_line_width(0.0)
        self.line(5.0, 5.0, 205.0, 5.0)  # top one
        self.line(5.0, 292.0, 205.0, 292.0)  # bottom one
        self.line(5.0, 5.0, 5.0, 292.0)  # left one
        self.line(205.0, 5.0, 205.0, 292.0)  # right one
    def add_pyplot (self,plt):
        temp_file_path = tempfile.NamedTemporaryFile(suffix=".png").name
        plt.savefig(temp_file_path, format='png', dpi=300)
        self.image(temp_file_path, w=150)
        os.remove(temp_file_path)
    def add_plotly(self,fig):
        temp_file_path = tempfile.NamedTemporaryFile(suffix=".png").name
        fig.write_image(temp_file_path)
        self.image(temp_file_path, w=150)
        os.remove(temp_file_path)
    def add_df (self, df):
        # Table Header
        self.add_page (orientation="landscape")
        ch =10
        self.set_text_color(0, 0, 0)
        self.set_font('Arial', 'B', 16)
        for c in df.columns:
            self.cell(40, h=ch, txt=c, border=1, align='C')

        # Table contents
        self.set_font('Arial', '', 16)
        for i, row in df.iterrows ():
            self.ln(ch)
            for c in df.columns:
                self.cell(w=40, h=ch,
                     txt=row[c].astype (str),
                     border=1,  align='C')

    def hyperlink(self, link_text, url):
        self.set_text_color(70, 150, 50)
        self.set_font("helvetica", "B", 16)
        self.cell (h=10, w=40, txt=link_text , link=url)

    def header(self):
        # Select Arial bold 15
        self.set_font('Arial', 'B', 15)
        if self.hdrimg != None:
            self.image (self.hdrimg, h=15, w = 15, y=10)
        if self.hdrtxt != None:
            self.cell(self.w, 15, self.hdrtxt, 0, 0, 'C')
        self.ln(20)
    def footer(self):
        # Go to 1.5 cm from bottom
        self.set_y(-15)
        # Select Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Print centered page number
        self.cell(0, 10, 'Page %s' % self.page_no(), 0, 0, 'C')

if __name__ == "__main__":
    xpdf = PDF_generator ('1859 Therapeutics','/Users/eric/OneDrive/Consulting/Projects/1859/euyb2djyortvcjnzxhpl.png' )
    xpdf.init_pdf()
    xpdf.border()

    xpdf.add_text ('Title', align='C', color=(255,0,0), fontsize = 40)
    xpdf.add_text('Next', bold='B')
    txt = 'There once was a man from Verdun\nNext line'

    xpdf.add_text(txt, bold='B')
    xpdf.hyperlink( "testlink", "https://cnn.com")
    xpdf.add_MoleculeFromSMILES (['CCCCCCCC', 'CCN'])
    import matplotlib.pyplot as plt

    data = {'a': np.arange(50),
            'c': np.random.randint(0, 50, 50),
            'd': np.random.randn(50)}
    data['b'] = data['a'] + 10 * np.random.randn(50)
    data['d'] = np.abs(data['d']) * 100

    plt.scatter('a', 'b', c='c', s='d', data=data)
    plt.xlabel('entry a')
    plt.xlim (0,50)
    plt.ylabel('entry b')
    plt.ylim(0, 50)
    xpdf.add_pyplot (plt)

    import plotly.graph_objects as go

    fig = go.Figure(data=
    go.Contour(
        z=[[10, 10.625, 12.5, 15.625, 20],
           [5.625, 6.25, 8.125, 11.25, 15.625],
           [2.5, 3.125, 5., 8.125, 12.5],
           [0.625, 1.25, 3.125, 6.25, 10.625],
           [0, 0.625, 2.5, 5.625, 10]],
        colorscale='Electric',
    ))
    xpdf.add_plotly(fig)
    df = pd.DataFrame ([[1,1,23], [3,5,67]],columns=['x','y','test'])
    xpdf.add_df (df)
    xpdf.save ()


#references
#https://towardsdatascience.com/creating-pdf-files-with-python-ad3ccadfae0f