#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 16:57:46 2022

@author: sookie
"""

import os
import argparse
#%%
parser = argparse.ArgumentParser()
parser.add_argument('-d','--dir', help='Path to directory containing the QTL csv files',required = True)
parser.add_argument('-l','--lod', help='LOD threshold',required = True)
parser.add_argument('-o','--out', help='Path to where you want to store the results',required = True)
parser.add_argument('-p1','--nam1', help='Parent1 name',required = True)
parser.add_argument('-p2','--nam2', help='Parent2 name',required = True)
parser.add_argument('-g1','--gff1', help='Path to the folder containing the GFF files',required = True)
parser.add_argument('-g2','--gff2', help='Path to the folder containing the GFF files',required = True)

args = parser.parse_args()
#%%
path = args.dir
out_dir = args.out
file_gff1 = args.gff1
file_gff2 = args.gff2
parent1 = args.nam2
parent2 = args.nam2
#%% ## ParentA

print("- Creation of", str(parent1)+" file with QTL regions + support interval + QTL trend")
with open(str(path)+str(parent1)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+"_with_trend.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\n")[0]
        line=line.split("\t")
#        print(line)
        chrm = line[0]
        binstart = str(line[1])
        binstop = str(line[2])
        peak_lod = float(line[3])
        peak_pos = line[4]
        interval = line[6]
        try:
            globals()[str(parent1)+'_support_interval_'+str(chrm)][str(binstart)+"-"+str(binstop)] = interval
        except:
            globals()[str(parent1)+'_support_interval_'+str(chrm)] = {}
#        break

#%%
print("- Creation of", str(parent1)+" file with QTL regions + support interval + QTL trend + gene overlap")
with open(str(path)+str(parent1)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+"_with_trend.csv","r") as scer, open(str(out_dir)+str(parent1)+"_QTL_regions_LOD_"+str(args.lod)+"_genes_overlap.csv","w") as out:
    out.write("Chrm\tBin Start (bp)\tBin Stop (bp)\tPeak LOD\tPeak (bp)\tLOD threshold\tLOD support interval\tQTL trend\tGene Name\tGene ID\tstart\tstop\tstrand\tWithin support interval?\n")
    for l in scer:
        l = l.split("\n")[0]
        l = l.split("\t")
        QTL_chrm = l[0]
        QTL_binstart = l[1]
        QTL_binstop = l[2]
        peak_lod = l[3]
        peak_bin = l[4]
        with open(str(file_gff1),"r") as file:
            for line in file:
                line = line.split("\n")[0]
                line = line.split("\t")
                chrm = line[0]
                nature = line[2]
                start = line[3]
                stop = line[4]
                strand = str(line[6])
                infos = line[8]
                geneID = infos.split(";")[0]
                geneName = infos.split(";")[1]
                if nature == "gene":
                    geneID = geneID.split("ID=")[1]
                    geneName = geneName.split("Name=")[1]
                    support = "no"
                    if chrm == QTL_chrm:
                        if float(start) >= float(QTL_binstart) and float(stop) <= float(QTL_binstop):
                            for k,v in globals()['scer_support_interval_'+str(chrm)].items():
                                qtl=k.split("-")
                                inter=v.split("-")
                                if float(qtl[0]) == float(QTL_binstart) and float(qtl[1]) == float(QTL_binstop):
                                    if float(start) >= float(inter[0]) and float(stop) <= float(inter[1]):
                                        support = "yes"
#                                        print("\nYES")
#                                        print("qtl,interval:",k,v)
#                                        print("qtl",binstart,binstop)
#                                        print("gene",start,stop)
                                    elif float(start) < float(inter[0]) and float(stop) < float(inter[0]):
#                                        print("\nNO")
#                                        print("qtl,interval:",k,v)
#                                        print("qtl",binstart,binstop)
#                                        print("gene",start,stop)
                                        support = "no"
                                    elif float(start) > float(inter[1]) and float(stop) > float(inter[1]):
#                                        print("\nNO")
#                                        print("qtl,interval:",k,v)
#                                        print("qtl",binstart,binstop)
#                                        print("gene",start,stop)
                                        support = "no"
                                    else:
                                        support = "partially overlap"
#                                        print("\nPartially")
#                                        print("qtl,interval:",k,v)
#                                        print("qtl",binstart,binstop)
#                                        print("gene",start,stop)
#                                    break
#                            print("gene",support,"overlap QTL",geneName,start,stop,"(QTL:",QTL_binstart,QTL_binstop,")")
                            out.write("\t".join(l)+"\t"+str(geneName)+"\t"+str(geneID)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(strand)+"\t"+str(support)+"\n")
                            
#%% ## skud
print("- Creation of", str(parent2)+" file with QTL regions + support interval + QTL trend")

with open(str(path)+str(parent2)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+"_with_trend.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\n")[0]
        line=line.split("\t")
#        print(line)
        chrm = line[0]
        binstart = str(line[1])
        binstop = str(line[2])
        peak_lod = float(line[3])
        peak_pos = line[4]
        interval = line[6]
        try:
            globals()[str(parent2)+'_support_interval_'+str(chrm)][str(binstart)+"-"+str(binstop)] = interval
        except:
            globals()[str(parent2)+'_support_interval_'+str(chrm)] = {}
#        break
#%%
print("- Creation of", str(parent2)+" file with QTL regions + support interval + QTL trend + gene overlap")
with open(str(path)+str(parent2)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+"_with_trend.csv","r") as skud, open(str(out_dir)+str(parent2)+"_QTL_regions_LOD_"+str(args.lod)+"_genes_overlap.csv","w") as out:
    out.write("Chrm\tBin Start (bp)\tBin Stop (bp)\tPeak LOD\tPeak (bp)\tLOD threshold\tLOD support interval\tQTL trend\tGene Name\tGene ID\tstart\tstop\tstrand\tWithin support interval?\n")
    for l in skud:
        l = l.split("\n")[0]
        l = l.split("\t")
        QTL_chrm = l[0]
        QTL_binstart = l[1]
        QTL_binstop = l[2]
        peak_lod = l[3]
        peak_bin = l[4]
        with open(str(file_gff2),"r") as file:
            for line in file:
                if line.startswith("#") == False:
                    line = line.split("\n")[0]
                    line = line.split("\t")
#                    print(line)
                    chrm = line[0]
                    nature = line[2]
                    start = line[3]
                    stop = line[4]
                    strand = str(line[6])
                    infos = line[8]
                    geneID = infos.split(";")[0]
                    geneName = infos.split(";")[1]
                    if nature == "gene":
                        geneID = geneID.split("ID=")[1]
                        geneName = geneName.split("Name=")[1]
                        support = "no"
                        if chrm == QTL_chrm:
                            if float(start) >= float(QTL_binstart) and float(stop) <= float(QTL_binstop):
                                for k,v in globals()[str(parent2)+'_support_interval_'+str(chrm)].items():
                                    qtl=k.split("-")
                                    inter=v.split("-")
                                    if float(qtl[0]) == float(QTL_binstart) and float(qtl[1]) == float(QTL_binstop):
                                        
                                        if float(start) >= float(inter[0]) and float(stop) <= float(inter[1]):
                                            support = "yes"
    #                                        print("\nYES")
    #                                        print("qtl,interval:",k,v)
    #                                        print("qtl",binstart,binstop)
    #                                        print("gene",start,stop)
                                        elif float(start) < float(inter[0]) and float(stop) < float(inter[0]):
    #                                        print("\nNO")
    #                                        print("qtl,interval:",k,v)
    #                                        print("qtl",binstart,binstop)
    #                                        print("gene",start,stop)
                                            support = "no"
                                        elif float(start) > float(inter[1]) and float(stop) > float(inter[1]):
    #                                        print("\nNO")
    #                                        print("qtl,interval:",k,v)
    #                                        print("qtl",binstart,binstop)
    #                                        print("gene",start,stop)
                                            support = "no"
                                        else:
                                            support = "partially overlap"
    #                                        print("\nPartially")
    #                                        print("qtl,interval:",k,v)
    #                                        print("qtl",binstart,binstop)
    #                                        print("gene",start,stop)
    #                                    break
    #                            print("gene",support,"overlap QTL",geneName,start,stop,"(QTL:",QTL_binstart,QTL_binstop,")")
                                out.write("\t".join(l)+"\t"+str(geneName)+"\t"+str(geneID)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(strand)+"\t"+str(support)+"\n")
                            