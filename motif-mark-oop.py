#!/usr/bin/env python
import cairo
import math
import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="specifies file inputs, index cutoff")
    parser.add_argument("-f", "--file", help="input: input FASTA file", \
        required=True, type=str)
    parser.add_argument("-i", "--int_fa", help="intermediate: one-line FASTA file", \
    required=True, type=str)
    parser.add_argument("-m", "--motifs", help="motif list file", \
        required=True, type=str)
    parser.add_argument("-o", "--outfile", help="output: png image of annotated genes file", \
        required=True, type=str)
    return parser.parse_args()

clinput = get_args()
f = clinput.file
i = clinput.int_fa
m = clinput.motifs
o = clinput.outfile

class Line():
    '''Line class contains subclasses'''
    def __init__(self, start, stop,r, kind, name=None, motif_kind=None):
        self.start=start
        self.stop=stop
        self.r=r
        self.kind=kind
        self.name=name
        self.motif_kind=motif_kind
        self.len=self.stop-self.start
        
        ## METHODS ##
    def __repr__(self):
        return f'{self.kind}({self.start}, {self.stop}, {self.motif_kind})'

    def draw(self):

        margin = 10
        width = 1000
        hgg = 100 #height of gene group
        height = hgg * r 

        
        if kind == "gene":
            context.set_line_width(3)
            context.set_source_rgb(0, 0, 0)
            context.move_to(margin+self.start, self.r*hgg)       
            context.line_to(margin+self.stop, self.r*hgg)
            context.set_font_size(10)
            context.move_to(margin+self.start, (self.r*hgg)-15) 
            context.show_text(self.name)
            context.stroke()
            surface.write_to_png(o)
        elif kind == "exon":
            context.set_line_width(8)
            context.set_source_rgb(255,0,0)
            context.move_to(margin+self.start, self.r*hgg)        #(x,y)
            context.line_to(margin+self.stop, self.r*hgg)
            context.stroke()
            surface.write_to_png(o)
        elif kind == "motif":
            context.set_line_width(10)
            context.set_source_rgb(red, green, blue)
            context.move_to(margin+self.start, self.r*hgg)        #(x,y)
            context.line_to(margin+self.stop, self.r*hgg)
            context.stroke()
            surface.write_to_png(o)
        
def oneline_fasta(FASTA, int_file):
    '''This function takes the sequence lines of a FASTA file and condenses them into a single line'''
    #!/usr/bin/env python
    seq = ''
    output = open(i, "w")
    with open(f, "r") as fh:
        for line in fh:
            line = line.strip('\n')
            if line.startswith('>'):
                if seq != "":
                    output.write(str(seq)+'\n')
                output.write(str(line)+'\n')
                seq=''
            else:
                seq+=line
        if seq:
            output.write(str(seq)+'\n')

def motifer(motif):
    '''this function transforms ambiguous bases into the bases they should be based on the IUPAC'''
    motif = motif.upper()
    amb_motif = motif.replace("U", "T")
    amb_motif = amb_motif.replace("W", "[AT]")
    amb_motif = amb_motif.replace("S", "[CG]")
    amb_motif = amb_motif.replace("M", "[AC]")
    amb_motif = amb_motif.replace("K", "[GT]")
    amb_motif = amb_motif.replace("R", "[AG]")
    amb_motif = amb_motif.replace("Y", "[CT]")
    amb_motif = amb_motif.replace("N", "[ACGT]")
    
    #look ahead
    # amb_motif = "(?=(" + amb_motif
    # amb_motif += "))"
    # motif = amb_motif
    return(amb_motif)


#Transforming input FASTA file into a one-line fasta file
oneline_fasta(f,i)

# Parsing Fasta File to create objects

#Initialize empty lists to contain ambiguous motifs and transformed motifs
motifs=[]
trans_motifs=[]

#color pallete for motifs, colors: green, yellow, blue, orange, fuschia
color_palette = [(0.4,0.7,0),(0.9,0.7,0), (0,0,1),(1,0,1),(0.6,0.2,0.8)]

#loop through motif file and save motifs into list, transform motifs and save them into another list
with open(m,'r') as mh:
    for line in mh:
        line=line.strip('\n')
        # print(line)
        motifs.append(line)
        trans_motif = motifer(line)
        trans_motifs.append(trans_motif)
        if line =='':
            break
    # print(motifs)
    # print(trans_motifs)

width, height = 1000, 500 
surface = cairo.PDFSurface("Figure_1.pdf",width, height)
context = cairo.Context(surface)

#Loop through one_line_fasta file, save genes, exons,and motifs as objects. Storing their start and stop positions, as well as motif sequence if object is a motif
#Save gene, exons, and motifs of a gene to a gene group so they can be drawn as a group
these_motifs=[]
r=0 #counter
with open(i,'r') as ih:
    for line in ih:
        line=line.strip('\n')
        if line.startswith('>'):
            header = line
            kind = "gene" #create string variable to define the class subtype
            r += 1
        else:
            gene = Line(0,len(line), r, kind, header) #save start, stop, type, and header as attributes in gene class
            gene.draw()
            start = 0 #initialize start as 0 so we can keep track of where the Upper case letter is (exon start)
            # print("starting seq enum of: ", header)
            for i, base in enumerate(line):
                if base.isupper():
                    if start != 0: #if start is not 0 then we are in an exon, continue until we encounter a lower case
                        continue
                    else:
                        start = i #holds position of the start of the exon if we have not encountered an uppercase yet
                elif base.islower():
                    if start != 0:
                        end = i-1 #once we encounter a lowercase, after starting an exon, save end position of the exon as the base preceeding the first lowercase
                        kind='exon' #subtype exon
                        exon = Line(start, end,r, kind)
                        exon.draw()
                        break
            # print(exon)            
            # print("Starting motif seeking:", len(motifs), len(line))
            for e, motif in enumerate(motifs):
                color = color_palette[e] #using enumerate, tie index of loop to color in color pallete list
                red=color[0] #set values in tupples of color_palette equal to variables so we can call them in draw
                green=color[1]
                blue=color[2]
                # print("searching for", motif)
                m = re.finditer(trans_motifs[e], line, re.IGNORECASE) #search for translated motifs indexed the same way as untranslated motifs
                # print("done")
                lc=0
                for k in m:
                    lc+=1
                    start, stop = (k.span())
                    # print(start, stop)
                    kind = 'motif'
                    this_motif = Line(start, stop, r, kind, motif_kind=motif)
                    this_motif.draw()
                    these_motifs.append(this_motif)
                    # print(len(these_motifs))
                # print("identified", lc, "motifs")


#Legend

#motif and height counters
mot = 925 #motif word x start position
start_hgt = 5 #legend item heigh start position
hgt =0 #height to increment by
col_start = 875 #color block start pos
col_stop = 885 #colior block end pos
context.rectangle(850,0,1000,150)
context.set_source_rgb(0, 0, 0) 
context.set_line_width(1)
context.stroke()   
for j,motif in enumerate(motifs):
    color = color_palette[j] #using enumerate, tie index of loop to color in color pallete list
    red=color[0] #set values in tupples of color_palette equal to variables so we can call them in draw
    green=color[1]
    blue=color[2]
    context.set_source_rgb(red, green, blue)
    hgt+=5
    context.set_line_width(10)
    context.move_to(col_start, start_hgt*hgt)       
    context.line_to(col_stop, start_hgt*hgt)
    
    context.set_font_size(10)
    context.move_to(mot,start_hgt*hgt)
    # print(hgt) 
    context.show_text(motif)
    context.stroke()

#Generating legend material for exons
context.set_source_rgb(1, 0, 0)
hgt+=5
context.set_line_width(10)
context.move_to(col_start, start_hgt*hgt)       
context.line_to(col_stop, start_hgt*hgt)
context.set_font_size(10)
context.move_to(mot,start_hgt*hgt) 
context.show_text("Exon")
context.stroke()
surface.write_to_png ("Figure_1.png")











