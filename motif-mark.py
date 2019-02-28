# Helena Klein
# Motif-marker, svg output, text and fasta input
# February 19, 2019
# Usage: python3 motif-mark.py -m <input motifs.txt> -g <input_genes.fasta> -r [if we want random colors(optional)]

import re
import math
import cairo
import random
import argparse

def get_arguments():
	parser = argparse.ArgumentParser(description="Produces an SVG figure of the gene region to scale with motifs highlighted. Gene of interest must have lowercase introns and uppercase exons. For updates, please contact helenasfklein@gmail.com")
	parser.add_argument('-m', '--motifs', required=True, help='Indicate the absolute path to a text file which contains all motifs of interest, one motif per line.')
	parser.add_argument('-g', '--genes', required=True, help='Indicate the absolute path to a fasta file which contains all genes of interest, needs to be 5 prime to 3 prime orientation.')
	parser.add_argument('-r', '--random', required=False, action='store_true', help='Generate random colors for motifs in figure')

	return parser.parse_args()
	
args = get_arguments()


# Argparse entries
motif_file=args.motifs
gene_file = args.genes

# For the color palette options (default is set colors)
rand=args.random

# Make sure the input is correct 
if ".fa" not in gene_file:
    exit("Gene file must be .fa or .fasta")
if ".txt" not in motif_file:
    exit("Motif file must be .txt")

# Get the full sequence out of the fasta per gene and make sure it's all in one string. 
# Save into a dictionary where the gene name is the key, return the dictionary

def get_genes(gene_file):
    genes = {}
    with open(gene_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if '>' in line:
                line = line.strip('>').split(' ')
                genes[line[0]] = ""
                gene = line[0]
            else:
                genes[gene] += line
            
    return genes

# Take the motifs of interest from a file and save in a set for easy reference. 
def get_motifs(motif_file):
    motifs = []
    with open(motif_file, "r") as fh:
        for line in fh:
            line = line.strip()
            motifs.append(line)
    motifs = set(motifs)
    return motifs

# Build a regex for a given motif (with ambiguity).
# These IUPAC identities came from this link: https://www.bioinformatics.org/sms/iupac.html

def build_regex(motif):
    regex = ""
    motif = motif.upper()
    regex_dict = {"A": "[Aa]", 
                 "T": "[Tt]",
                 "C": "[Cc]",
                 "G": "[Gg]",
                 "N": "[AaTtGgCc]",
                 "Y": '[ctCT]',
                 "R": '[AaGg]',
                 "S": '[GgCc]',
                 "W": '[AaTt]',
                 "K": '[GgTt]',
                 "M": '[AaCc]',
                 "B": '[CcGgTt]',
                 "D": '[AaGgTt]',
                 "H": "[AaCcTt]",
                 "V": "[AaCcGg]"}
    
    for position in range(len(motif)):
        base = motif[position]
        if base in regex_dict:
            regex += regex_dict[base]
        else:
            exit("Invalid character input in motif:", motif)
    return regex


# Here is the main body of my code now that I have defined my functions

# Take input files and convert them to useable formats
genes = get_genes(gene_file)
motifs = get_motifs(motif_file)
if len(motifs) > 8:
    exit("Too many motifs, maximum is 8")


# Break each gene into kmers for each motif and record where the motif occurred in these genes
gene_motif = {}
motif_loc = {}
exon_loc={}
for gene in genes:
    # Make sure to record on a one based counting system
    sequence = genes[gene]
    
    for motif in motifs:
        # record positions of found motifs
        
        length = len(motif)
        # get the regex for the motif as needed
        regex = build_regex(motif)
        
        motif_loc[motif] = []
        
        for position in range(len(sequence)):
            kmer =  sequence[int(position):position+length]
            if re.match(regex, kmer):
                position += 1
                motif_loc[motif].append(position)
    gene_motif[gene] = motif_loc
        
        # For exon locations, if the position isn't one away from the last stored position,
        #then record the previous (end of last exon) and the current position (start of next 
        #exon) and save the last one as the end of them when we reach the end of the gene
        
    
    exon_number = 0
    first=True
    start = 0
    counter = 0
    exon = False
    
      
    exon_loc[gene] = {}    
    for position in range(len(sequence)):
        if sequence[position].isupper():
            position += 1
            
            # record start of exon
            if first:
                exon_number += 1
                
                exon_loc[gene][exon_number] =[]
                start = position 
                first=False
                exon=True
            counter += 1
        else:
            first = True
            # Stick locations of each exon for each gene into a dictionary of dictionaries
            if exon:
                exon_loc[gene][exon_number].append(start)
                exon_loc[gene][exon_number].append(start+counter-1)
                
                counter=0
                exon=False
            
                    
# Scale the locations so they are easy to use in image drawing
image_stats = {}
for gene in genes:
    length = len(genes[gene])
    scale = (1/length) * .9
    shift = 0.05
    
    # Scale motif locations
    for motif in gene_motif[gene]:
        loc = []
        for location in gene_motif[gene][motif]:
            loc.append((location*scale)+shift)
        gene_motif[gene][motif] = loc
    # Scale exon locations
    for exon in exon_loc[gene]:
        loc = []
        for location in exon_loc[gene][exon]:
            loc.append((location*scale)+shift)
        exon_loc[gene][exon] = loc
    


# Make an image using pycairo
width, height = 800, 500
name = "{}.svg"

for gene in genes:
    #create the coordinates to display your graphic, designate output
    file = name.format(gene)
    surface = cairo.SVGSurface(file,width, height)
    #create the coordinates you will be drawing on (like a transparency) - you can create a transformation matrix
    context = cairo.Context(surface)
    context.scale(width,height) #will set your drawing surface to a 0.0-1.0 scale
    
    # Add a title
    context.select_font_face("Times", cairo.FONT_SLANT_NORMAL, 
    cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(.1)
    context.move_to(shift, 0.25)
    context.show_text(gene)
    
    # Figure Legend
    context.set_font_size(.06)
    context.move_to(0.8, 0.25)
    context.show_text("Motifs")
    
    # exon legend entry
    context.select_font_face("Calibri", cairo.FONT_SLANT_NORMAL, 
    cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(.03)
    context.move_to(0.8, 0.30)
    legend_pos = 0.30
    context.show_text("Exon")
    
    #Need to tell cairo where to put the brush, the color and width, and the shape you want it to draw
    #draw a line
    context.set_line_width(0.005)
    context.move_to(shift,0.75)        #(x,y)
    context.line_to(1-shift,0.75)
    context.stroke()
    
    # Draw on exons
    for exon in exon_loc[gene]:
        start = exon_loc[gene][exon][0]
        width = exon_loc[gene][exon][1]-start
        
        # rectangle
        context.rectangle(start,0.725,width,0.05)
        context.fill()
    motif_number = 0
    palette =  (.8,0.1,0.1), (0.1,0.9,0.9), (.4, 0.9, 0.2), (.6,0.2,0.8), (.0, 0.0, 0.9), (.9, 0.9,0), (0.9,0.3,0.8), (.9,.6,0)
    for motif in gene_motif[gene]:
        # Make color palette
        if rand:
            r = round(random.uniform(0,1),1)
            g = round(random.uniform(0,1),1)
            b = round(random.uniform(0,1),1)
            
            context.set_source_rgb(r, g, b)
        else:
            r = palette[motif_number][0]
            g = palette[motif_number][1]
            b = palette[motif_number][2]
            
            context.set_source_rgb(r,g,b)
            motif_number += 1
        
        entry = motif.upper()
        
        legend_pos += 0.05
        
        # legend entry
        context.set_font_size(.03)
        context.move_to(0.8, legend_pos)
        context.show_text(entry)
        
        # rectangle
        for location in gene_motif[gene][motif]:
            
            context.rectangle(location,0.725,0.002,0.05)
            context.fill()
        
    
    
    #write out drawing
    surface.finish()



