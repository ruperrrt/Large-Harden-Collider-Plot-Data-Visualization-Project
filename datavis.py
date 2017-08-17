#!/usr/bin/env python
import math
import cairocffi as cairo
img = {}  
img['height'] = 1000
img['width'] = 1000
img['center_x']   = img['width']  / 2.0
img['center_y']   = img['height'] / 2.0
img['font_size']  = 16
# store fst and rna 
in_file1="/home/a-m/ib501_stud24/shell/hw11/Gacu_FoldChange_GenomCoords.tsv"
fh_1=open(in_file1, "r")
rna_stats = {}
# Add each chromosome to the dictionary and store the
# basepair and statistical value
for line in fh_1:
    lines=line.split()
    if lines[1].startswith('group'):
        if lines[1] in rna_stats: 
            rna_stats[lines[1]]['(bp,stat)'].append((int(lines[2]),float(lines[4])))
        if lines[1] not in rna_stats:
            rna_stats[lines[1]]={}
            rna_stats[lines[1]]['(bp,stat)']=[]
            rna_stats[lines[1]]['(bp,stat)'].append((int(lines[2]),float(lines[4])))
        rna_stats[lines[1]]['(bp,stat)'].sort()
#
#
# I use tuple here because it is comparable and can be sorted by bp
#

in_file2="/home/a-m/ib501_stud24/shell/hw11/batch_1.phistats_fw2-oc.tsv"
fh_2=open(in_file2, "r")
fst_stats = {}
# Add each chromosome to the dictionary and store the
# basepair and statistical value
for line in fh_2:
    lines=line.split()
    if len(lines)>4:
        if lines[4].startswith('group'):
            if lines[4] in fst_stats: 
                fst_stats[lines[4]]['(bp,stat)'].append((int(lines[5]),float(lines[10])))
            if lines[4] not in fst_stats:
                fst_stats[lines[4]]={}
                fst_stats[lines[4]]['(bp,stat)']=[]
                fst_stats[lines[4]]['(bp,stat)'].append((int(lines[5]),float(lines[10])))
                fst_stats[lines[4]]['(bp,stat)'].sort()

 
ps = cairo.PDFSurface('/home/a-m/ib501_stud24/shell/hw11/Outfile' , img['height'], img['width'])
cr = cairo.Context(ps)
#all compution starts from here
# Convert a radius and a span of degrees into X, Y coordinates
def get_x_y_coordinates(center_x, center_y, degree, radius):
    if degree <= 90:
        theta=float(degree)
        opp_side=radius * math.sin(math.radians(theta))
        adj_side=radius * math.cos(math.radians(theta))
        x=center_x + adj_side
        y=center_y + opp_side
    elif degree <= 180:
        theta=float(degree-90.0)
        opp_side=radius * math.sin(math.radians(theta))
        adj_side=radius * math.cos(math.radians(theta))
        x=center_x - opp_side
        y=center_y + adj_side
    elif degree <= 270:
        theta=float(degree-180.0)
        opp_side=radius * math.sin(math.radians(theta))
        adj_side=radius * math.cos(math.radians(theta))
        x=center_x - adj_side
        y=center_y - opp_side
    else:
        theta=float(degree-270.0)
        opp_side=radius * math.sin(math.radians(theta))
        adj_side=radius * math.cos(math.radians(theta))
        x=center_x + opp_side
        y=center_y - adj_side
    return x,y

# store the dimension  for each chromosome in a dictionary
chrom_size={'groupI'    : 28185914, 'groupII'   : 23295652, 'groupIII'   : 16798506, 'groupIV'  : 32632948,
            'groupV'    : 12251397, 'groupVI'    : 17083675, 'groupVII' : 27937443,
            'groupVIII' : 19368704, 'groupIX'   : 20249479,'groupX'    : 15657440, 'groupXI'    : 16706052, 'groupXII' : 18401067,
            'groupXIII' : 20083130, 'groupXIV'  : 15246461  , 'groupXV'  : 16198764,
            'groupXVI'  : 18115788, 'groupXVII' : 14603141, 'groupXVIII' : 16282716, 'groupXIX'   : 20240660,'groupXX'  : 19732071,
            'groupXXI'  : 11717487}
key=['groupI','groupII','groupIII','groupIV','groupV','groupVI','groupVII','groupVIII','groupIX','groupX','groupXI','groupXII','groupXIII','groupXIV','groupXV','groupXVI','groupXVII','groupXVIII','groupXIX','groupXX','groupXXI']


 # create a list to store chromosome number in an increasing order
total=0
for keys in key:
    total+=chrom_size[keys]
for keys in key:
    chrom_size[keys]=chrom_size[keys]*1.0/total  # to update value with percentage
totaldegree=360-len(chrom_size)*2 # set gap in graph as 2 degree 

#create a function to set hsv color for inner circle
def color_in(value):
    if value>=0.9 and value<=1:
        x,y,z=0.66,0.21,0.21
    if value>=0.8 and value<0.9:
        x,y,z=0.66,0.29,0.21
    if value>=0.7 and value<0.8:
        x,y,z=0.66,0.37,0.21
    if value>=0.6 and value<0.7:
        x,y,z=0.66,0.44,0.21
    if value>=0.5 and value<0.6:
        x,y,z=0.66,0.51,0.21
    if value>=0.4 and value<0.5:
        x,y,z=0.66,0.59,0.21
    if value>=0.3 and value<0.4:
        x,y,z=0.66,0.66,0.21
    if value>=0.2 and value<0.3:
        x,y,z=0.59,0.66,0.21
    if value>=0.1 and value<0.2:
        x,y,z=0.51,0.66,0.21
    if value<0.1:
        x,y,z=0.44,0.66,0.21
    return x,y,z

#create a function to set hsv color for outter circle
def color_out(value):
    if value>=0.9 and value<=1:
        x,y,z=0.55,0,0.66
    if value>=0.8 and value<0.9:
        x,y,z=0.44,0,0.66
    if value>=0.7 and value<0.8:
        x,y,z=0.33,0,0.66
    if value>=0.6 and value<0.7:
        x,y,z=0.22,0,0.66
    if value>=0.5 and value<0.6:
        x,y,z=0.11,0,0.66
    if value>=0.4 and value<0.5:
        x,y,z=0,0,0.66
    if value>=0.3 and value<0.4:
        x,y,z=0,0.11,0.66
    if value>=0.2 and value<0.3:
        x,y,z=0,0.22,0.66
    if value>=0.1 and value<0.2:
        x,y,z=0,0.33,0.66
    if value<0.1:
        x,y,z=0,0.44,0.66
    return x,y,z

# drawing the key for div starts from here
ind=[0,1,2,3,4,5,6,7,8,9]
cr.set_line_width(0.01)
for val in ind:
    cr.rectangle(800,30+val*15,40,15)
    x,y,z=color_out(1-0.1*(val+1))
    cr.set_source_rgb(x, y, z)
    cr.fill()
cr.set_line_width(1)
cr.set_source_rgb(0, 0, 0)
cr.rectangle(800, 30, 40, 150)
cr.stroke()
textents = cr.text_extents('Div')
text_width  = textents[2]
text_height = textents[3]
x,y =820,30
cr.move_to(x-text_width/2.0, y-text_height/2.0)
cr.show_text('Div')
textents = cr.text_extents('1.0')
text_width  = textents[2]
text_height = textents[3]
x,y =790,40
cr.move_to(x-text_width/2.0, y-text_height/2.0)
cr.show_text('1.0')
text_width  = textents[2]
text_height = textents[3]
x,y =790,115
cr.move_to(x-text_width/2.0, y-text_height/2.0)
cr.show_text('0.5')
text_width  = textents[2]
text_height = textents[3]
x,y =790,180
cr.move_to(x-text_width/2.0, y-text_height/2.0)
cr.show_text('0')

# drawing the key for fst starts from here
ind=[0,1,2,3,4,5,6,7,8,9]
cr.set_line_width(0.01)
for val in ind:
    cr.rectangle(900,30+val*15,40,15)
    x,y,z=color_in(1-0.1*(val+1))
    cr.set_source_rgb(x, y, z)
    cr.fill()
cr.set_line_width(1)
cr.set_source_rgb(0, 0, 0)
cr.rectangle(900, 30, 40, 150)
cr.stroke()
textents = cr.text_extents('Fst')
text_width  = textents[2]
text_height = textents[3]
x,y =920,30
cr.move_to(x-text_width/2.0, y-text_height/2.0)
cr.show_text('Fst')
textents = cr.text_extents('1.0')
text_width  = textents[2]
text_height = textents[3]
x,y =890,40
cr.move_to(x-text_width/2.0, y-text_height/2.0)
cr.show_text('1.0')
text_width  = textents[2]
text_height = textents[3]
x,y =890,115
cr.move_to(x-text_width/2.0, y-text_height/2.0)
cr.show_text('0.5')
text_width  = textents[2]
text_height = textents[3]
x,y =890,180
cr.move_to(x-text_width/2.0, y-text_height/2.0)
cr.show_text('0')

# drawing the dashed lines starts from here
startdegree=0
cr.set_dash([6.0,1.5])
cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
cr.set_font_size(8)
textents = cr.text_extents('10Mb')
text_width  = textents[2]
text_height = textents[1]
for keys in key:
    num=chrom_size[keys]*total/5000000   #calculate the number of dashed lines per chromosome
    num=int(num)
    i=1
    degree=startdegree
    while(num>0):
        degree=degree+1.0*totaldegree*5000000/total
        x1,y1=get_x_y_coordinates(img['center_x'], img['center_y'],degree,250)
        x2,y2=get_x_y_coordinates(img['center_x'], img['center_y'],degree,350)
        x3,y3=get_x_y_coordinates(img['center_x'], img['center_y'],degree,235)
        cr.set_source_rgb(0.5,0.5,0.5)
        if i%2 ==1 :
            cr.set_line_width(1.5)
            cr.move_to(x1,y1)
            cr.line_to(x2,y2)
            cr.stroke()
        else:
            cr.set_line_width(4)
            cr.move_to(x1,y1)
            cr.line_to(x2,y2)
            cr.stroke()
            if i == 2:
                cr.set_source_rgb(0, 0, 0)
                textents = cr.text_extents('10Mb')
                cr.move_to(x3-text_width/2.0, y3-text_height/2.0)
                cr.show_text('10Mb')
            if i == 4:
                cr.set_source_rgb(0, 0, 0)
                textents = cr.text_extents('20Mb')
                cr.move_to(x3-text_width/2.0, y3-text_height/2.0)
                cr.show_text('20Mb')
            if i == 6:
                cr.set_source_rgb(0, 0, 0)
                textents = cr.text_extents('30Mb')
                cr.move_to(x3-text_width/2.0, y3-text_height/2.0)
                cr.show_text('30Mb')
        num=num-1
        i=i+1
    startdegree=startdegree+2+totaldegree*chrom_size[keys]
cr.move_to(770,500)
cr.set_dash([1.0,0.0])
cr.set_line_width(2)
# drawing the inner circle starts from here
i=0        # i is the starting degree of each arc
for keys in key: 
    cr.arc(img['center_x'], img['center_y'], 270, math.radians(i), math.radians(i+totaldegree*chrom_size[keys]))
    cr.arc_negative(img['center_x'], img['center_y'], 300, math.radians(i+totaldegree*chrom_size[keys]), math.radians(i))
    cr.close_path()
    cr.set_source_rgb(0, 0, 0)
    cr.stroke_preserve()
    cr.set_source_rgb(0.5, 0.5, 0.5)
    cr.fill()
    i=i+2+totaldegree*chrom_size[keys]   # set the starting degree of next arc

# drawing the outter circle starts from here

i=0        # i is the starting degree of each arc
for keys in key:
    cr.arc(img['center_x'], img['center_y'], 310, math.radians(i), math.radians(i+totaldegree*chrom_size[keys]))
    cr.arc_negative(img['center_x'], img['center_y'], 340, math.radians(i+totaldegree*chrom_size[keys]), math.radians(i))
    cr.close_path()
    cr.set_source_rgb(0, 0, 0)
    cr.stroke_preserve()
    cr.set_source_rgb(0.5, 0.5, 0.5)
    cr.fill()
    i=i+2+totaldegree*chrom_size[keys]   # set the starting degree of next arc

#coloring the inner circle starts from here

i=0        # i is the starting degree of each arc
for keys in key:
    j=0
    start=i
    while j<len(fst_stats[keys]['(bp,stat)'])-1:
        cr.arc(img['center_x'], img['center_y'], 270, math.radians(start), math.radians(start+totaldegree*(fst_stats[keys]['(bp,stat)'][j+1][0]-fst_stats[keys]['(bp,stat)'][j][0])*1.0/total))
        cr.arc_negative(img['center_x'], img['center_y'], 300, math.radians(start+totaldegree*(fst_stats[keys]['(bp,stat)'][j+1][0]-fst_stats[keys]['(bp,stat)'][j][0])*1.0/total), math.radians(start))
        cr.close_path()
        x,y,z=color_in(fst_stats[keys]['(bp,stat)'][j][1])
        cr.set_source_rgb(x,y,z)
        cr.stroke_preserve()
        cr.set_source_rgb(x,y,z)
        cr.fill()
        start=start+totaldegree*(fst_stats[keys]['(bp,stat)'][j+1][0]-fst_stats[keys]['(bp,stat)'][j][0])*1.0/total
        j=j+1
    i=i+2+totaldegree*chrom_size[keys]

#coloring the outter circle starts from here

i=0        # i is the starting degree of each arc
for keys in key:
    j=0
    start=i
    while j<len(rna_stats[keys]['(bp,stat)'])-1:
        cr.arc(img['center_x'], img['center_y'], 310, math.radians(start), math.radians(start+totaldegree*(rna_stats[keys]['(bp,stat)'][j+1][0]-rna_stats[keys]['(bp,stat)'][j][0])*1.0/total))
        cr.arc_negative(img['center_x'], img['center_y'], 340, math.radians(start+totaldegree*(rna_stats[keys]['(bp,stat)'][j+1][0]-rna_stats[keys]['(bp,stat)'][j][0])*1.0/total), math.radians(start))
        cr.close_path()
        x,y,z=color_out(rna_stats[keys]['(bp,stat)'][j][1])
        cr.set_source_rgb(x,y,z)
        cr.stroke_preserve()
        cr.set_source_rgb(x,y,z)
        cr.fill()
        start=start+totaldegree*(rna_stats[keys]['(bp,stat)'][j+1][0]-rna_stats[keys]['(bp,stat)'][j][0])*1.0/total
        j=j+1
    i=i+2+totaldegree*chrom_size[keys]

# adding txt for outter circle
cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)
cr.set_font_size(12)
cr.set_source_rgb(0, 0, 0)
# textents = cr.text_extents(keys[5:])
startdegree=0
for keys in key:
    textents = cr.text_extents(keys[5:])                                                                               
    text_width  = textents[2]
    text_height = textents[3]
    x,y = get_x_y_coordinates(img['center_x'], img['center_y'], startdegree+totaldegree*chrom_size[keys]/2 , 380)
    cr.move_to(x-text_width/2.0, y-text_height/2.0)
    cr.show_text(keys[5:])
    startdegree=startdegree+2+totaldegree*chrom_size[keys]
#cr.set_line_width(2)
# cr.close_path()
# Close the file
cr.show_page()
