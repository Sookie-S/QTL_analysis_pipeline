#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 14:39:09 2022

@author: sookie
"""

#%%
import os
import argparse
#%%
parser = argparse.ArgumentParser()
parser.add_argument('-l','--lod', help='',required = True)
parser.add_argument('-i','--int', help='',required = True)
parser.add_argument('-d','--dir', help='Path to the Input directory containing the Files',required = True)
parser.add_argument('-o','--out', help='Path to the Output directory',required = True)
parser.add_argument('-p1','--pre1', help='Chromosome prefix parent 1',required = True)
parser.add_argument('-p2','--pre2', help='Chromosome prefix parent 2',required = True)
parser.add_argument('-n1','--nam1', help='Name parent 1',required = True)
parser.add_argument('-n2','--nam2', help='Name parent 2',required = True)

args = parser.parse_args()
#%%
path = args.dir
lod_threshold = args.lod
lod_threshold = float(lod_threshold)
interval = args.int
interval = float(interval)
out_dir = args.out
prefix1 = args.pre1
prefix2 = args.pre2
name1 = args.nam1
name2 = args.nam2
#%%
folders = os.listdir(path)
for f in folders:
    if str(name2) in str(f):
        files_parent2 = os.listdir(str(path)+"/"+str(f))
        folders_parent2 = f
    elif str(name1) in str(f):
        files_parent1 = os.listdir(str(path)+"/"+str(f))
        folders_parent1 = f
#%%
with open(str(out_dir)+"tmp_file_"+str(name1)+"_QTL_regions_details_LOD_"+str(args.lod)+".csv","w") as out, open(str(out_dir)+str(name1)+"_QTL_regions_LOD_"+str(args.lod)+".csv","w") as out2:
    out.write("Chrm\tBin start (bp)\tMLE allele freq.\tLOD score\n")
    out2.write("Chrm\tBin Start (bp)\tBin Stop (bp)\tPeak LOD\tPeak (bp)\tLOD threshold\n")
    for f in files_parent1:
        if ".out" in f:
            chrm = f.split(".out")[0]
            chrm = chrm.split("_")[1]
            print("- Current file:",f,chrm)
            globals()[str(name1)+'_lod_dico_'+str(chrm)] = {}
            with open(str(path)+str(folders_parent1)+"/"+str(f),"r") as file:
                next(file)
                for line in file:
                    ok = False
                    line = line.split("\n")[0]
                    line = line.split("\t")
                    binstart = line[0]
                    lod = line[2]
                    globals()[str(name1)+'_lod_dico_'+str(chrm)][binstart] = lod
                    if float(lod) >= float(lod_threshold):
                        ok = True
                        bin_start = binstart
                        bin_stop = binstart
                        peak = binstart
                        peak_lod = lod
#                        print("\nBin start:",bin_start,lod,chrm)
                        out.write(str(chrm)+"\t"+"\t".join(line)+"\n")
                        while ok:
                            try:
                                line = next(file)
                            except:
#                                print("END OF FILE")
                                break
                            line = line.split("\n")[0]
                            line = line.split("\t")
                            binstart = line[0]
                            lod = line[2]
                            globals()[str(name1)+'_lod_dico_'+str(chrm)][binstart] = lod
                            if float(lod) >= float(lod_threshold):
                                ok = True
                                out.write(str(chrm)+"\t"+"\t".join(line)+"\n")
                                bin_stop = binstart
                                if float(lod) > float(peak_lod):
                                    peak_lod = lod
                                    peak = binstart
                            else:
#                                print("Bin stop:",bin_stop,lod,chrm)
                                if float(bin_stop)-float(bin_start) >= 20000:
                                    out2.write(str(chrm)+"\t"+str(bin_start)+"\t"+str(bin_stop)+"\t"+str(peak_lod)+"\t"+str(peak)+"\t"+str(lod_threshold)+"\n")
                                ok = False
                                out.write(" \t \t \t \n")
                        else:
#                            print("end WHILE")
#                            print(line)
                            ok = False
#%%
with open(str(out_dir)+"tmp_file_"+str(name2)+"_QTL_regions_details_LOD_"+str(args.lod)+".csv","w") as out, open(str(out_dir)+str(name2)+"_QTL_regions_LOD_"+str(args.lod)+".csv","w") as out2:
    out.write("Chrm\tBin start (bp)\tMLE allele freq.\tLOD score\n")
    out2.write("Chrm\tBin Start (bp)\tBin Stop (bp)\tPeak LOD\tPeak (bp)\tLOD threshold\n")
    for f in files_parent2:
        if ".out" in f:
            chrm = f.split(".out")[0]
            chrm = chrm.split("_",1)[1]
            print("- Current file:",f,chrm)
            globals()[str(name2)+'_lod_dico_'+str(chrm)] = {}
            with open(str(path)+str(folders_parent2)+"/"+str(f),"r") as file:
                next(file)
                for line in file:
                    ok = False
#                    print("new",line)
                    line = line.split("\n")[0]
                    line = line.split("\t")
                    binstart = line[0]
                    lod = line[2]
                    globals()[str(name2)+'_lod_dico_'+str(chrm)][binstart] = lod
                    if float(lod) >= float(lod_threshold):
                        ok = True
                        bin_start = binstart
                        bin_stop = binstart
                        peak = binstart
                        peak_lod = lod
#                        print("\nBin start:",bin_start,lod,chrm)
                        out.write(str(chrm)+"\t"+"\t".join(line)+"\n")
                        while ok:
#                            print("OK2")
#                            print(line)
                            try:
                                line = next(file)
                            except:
#                                print("End of file")
                                break
                            line = line.split("\n")[0]
                            line = line.split("\t")
                            binstart = line[0]
                            lod = line[2]
                            globals()[str(name2)+'_lod_dico_'+str(chrm)][binstart] = lod
                            if float(lod) >= float(lod_threshold):
#                                print("YEES")
                                ok = True
                                out.write(str(chrm)+"\t"+"\t".join(line)+"\n")
                                bin_stop = binstart
                                if float(lod) > float(peak_lod):
                                    peak_lod = lod
                                    peak = binstart
                            else:
#                                print("Bin stop:",bin_stop,lod,chrm)
                                if float(bin_stop)-float(bin_start) >= 20000:
                                    out2.write(str(chrm)+"\t"+str(bin_start)+"\t"+str(bin_stop)+"\t"+str(peak_lod)+"\t"+str(peak)+"\t"+str(lod_threshold)+"\n")
                                ok = False
                                out.write(" \t \t \t \n")
                        else:
#                            print("end WHILE")
#                            print(line)
                            ok = False

#%% add lod support interval
#%%
with open(str(out_dir)+str(name1)+"_QTL_regions_LOD_"+str(args.lod)+".csv","r") as file, open(str(out_dir)+str(name1)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+".csv","w") as out:
    header = next(file)
    header = header.split("\n")[0]
    header = header+"\tLOD support interval\n"
    out.write(header)
    for line in file:
        line = line.split("\n")[0]
        line = line.split("\t")
        chrm = line[0]
        binstart = float(line[1])
        binstop = float(line[2])
        peak_lod = float(line[3])
        peak_pos = line[4]
#        print(line)

#        print("Peak:",peak_lod,"at the position",peak_pos)
        peak_min = peak_lod - interval
        peak_max = peak_lod + interval
#        print("Lod interval:",peak_min,peak_max)
        support_int_start = ""
        support_int_stop = ""
        ok = False
        for k,v in globals()[str(name1)+'_lod_dico_'+str(chrm)].items():
            if float(k) >= binstart and float(k) <= binstop:
                if float(v) >= peak_min and float(v) <= peak_max:
                    if ok == False:
                        support_int_start = k
                        ok = True
                    support_int_stop = k
    #                print("within peak interval",k,v)
#        print("QTL",binstart,binstop)
#        print("support interval",support_int_start,support_int_stop)
        out.write("\t".join(line)+"\t"+str(support_int_start)+"-"+str(support_int_stop)+"\n")
#        break
#%%
with open(str(out_dir)+str(name2)+"_QTL_regions_LOD_"+str(args.lod)+".csv","r") as file, open(str(out_dir)+str(name2)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+".csv","w") as out:
    header = next(file)
    header = header.split("\n")[0]
    header = header+"\tLOD support interval\n"
    out.write(header)
    for line in file:
        line = line.split("\n")[0]
        line = line.split("\t")
        chrm = line[0]
        binstart = float(line[1])
        binstop = float(line[2])
        peak_lod = float(line[3])
        peak_pos = line[4]
#        print(line)

#        print("Peak:",peak_lod,"at the position",peak_pos)
        peak_min = peak_lod - interval
        peak_max = peak_lod + interval
#        print("Lod interval:",peak_min,peak_max)
        support_int_start = ""
        support_int_stop = ""
        ok = False
        for k,v in globals()[str(name2)+'_lod_dico_'+str(chrm)].items():
            if float(k) >= binstart and float(k) <= binstop:
                if float(v) >= peak_min and float(v) <= peak_max:
                    if ok == False:
                        support_int_start = k
                        ok = True
                    support_int_stop = k
    #                print("within peak interval",k,v)
#        print("QTL",binstart,binstop)
#        print("support interval",support_int_start,support_int_stop)
        out.write("\t".join(line)+"\t"+str(support_int_start)+"-"+str(support_int_stop)+"\n")
#        break

#%%

## search trend QTL

#%%

#%% ## parent 1
with open(str(out_dir)+"tmp_file_"+str(name1)+"_QTL_regions_details_LOD_"+str(args.lod)+".csv","r") as file:
    next(file)
    line=next(file)
    line=line.split("\n")[0]
    line=line.split("\t")
    chrm = str(line[0])
    binstart = str(line[1])
    alleleFreq = str(line[2])
    LOD = line[3]
    # print(line)
#    print("\n--START OF QTL REGION--")
    trend = {}
    status = "Up"
    QTLtrend = []
    QTLtrend.append(status)
    for line in file:
        if str(prefix1) not in chrm:
#            print("\n--START OF QTL REGION--")
            trend[str(previous_chrm)+"_"+str(previous_binstart)]=QTLtrend
            QTLtrend = []
            status = "Up"
            QTLtrend.append(status)
            chrm = str(prefix1)
            binstart = 0
            alleleFreq = 0
            LOD = 0
        # print(line)
        previousLOD = LOD
        previous_chrm = chrm
        previous_binstart = binstart
        line=line.split("\n")[0]
        line=line.split("\t")
        chrm = str(line[0])
        binstart = str(line[1])
        alleleFreq = str(line[2])
        LOD = line[3]
        if str(prefix1) in chrm and str(prefix1) in previous_chrm:
            # print(line)
            # print("values:",LOD,previousLOD)
            if float(LOD) > float(previousLOD):
#                print("Up")
                newstatus = "Up"
            elif float(LOD) < float(previousLOD):
#                print("Down")
                newstatus = "Down"
            else:
#                print("equal")
                newstatus = status
            
        if status != newstatus:
            status = newstatus
            QTLtrend.append(status)
            
    try:
        next(line)
    except:
        trend[str(previous_chrm)+"_"+str(previous_binstart)]=QTLtrend

#%%
with open(str(out_dir)+str(name1)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+".csv","r") as file, open(str(out_dir)+str(name1)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+"_with_trend.csv","w") as out:
    out.write("Chrm\tBin Start (bp)\tBin Stop (bp)\tPeak LOD\tPeak (bp)\tLOD threshold\tLOD support interval\tQTL trend\n")
    next(file)
    for line in file:
        ok = False
        line=line.split("\n")[0]
        line=line.split("\t")
#        print("\n",line)
        chrm = line[0]
        binstop = line[2]
#        print(binstop)
        for k,v in trend.items():
            k_chr = k.split("_")[0]
            k_binstop = k.split("_")[1]
            if k_binstop == binstop and k_chr == chrm:
                ok = True
#                print("  --YEEES",k,v)
                out.write("\t".join(line)+"\t"+"-".join(v)+"\n")
        if ok == False:
            print("  --NOT FOUND")
            break
#%%
#%% ## parent 2
with open(str(out_dir)+"tmp_file_"+str(name2)+"_QTL_regions_details_LOD_"+str(args.lod)+".csv") as file:
    next(file)
    line=next(file)
    line=line.split("\n")[0]
    line=line.split("\t")
    chrm = str(line[0])
    binstart = str(line[1])
    alleleFreq = str(line[2])
    LOD = line[3]
#    print(line)
#    print("\n--START OF QTL REGION--")
    trend = {}
    status = "Up"
    QTLtrend = []
    QTLtrend.append(status)
    for line in file:
        if str(prefix2) not in chrm:
#            print(str(previous_chrm)+"-"+str(previous_binstart))
#            print(QTLtrend)
#            print("\n--START OF QTL REGION--")
            trend[str(previous_chrm)+"-"+str(previous_binstart)]=QTLtrend
            QTLtrend = []
            status = "Up"
            QTLtrend.append(status)
            chrm = str(prefix2)
            binstart = 0
            alleleFreq = 0
            LOD = 0
#        print(line)
        previousLOD = LOD
        previous_chrm = chrm
        previous_binstart = binstart
        line=line.split("\n")[0]
        line=line.split("\t")
        chrm = str(line[0])
        binstart = str(line[1])
        alleleFreq = str(line[2])
        LOD = line[3]
        if str(prefix2) in chrm and str(prefix2) in previous_chrm:
            # print(line)
            # print("values:",LOD,previousLOD)
            if float(LOD) > float(previousLOD):
#                print("Up")
                newstatus = "Up"
            elif float(LOD) < float(previousLOD):
#                print("Down")
                newstatus = "Down"
            else:
#                print("equal")
                newstatus = status
            
        if status != newstatus:
            status = newstatus
            QTLtrend.append(status)
    try:
        next(line)
    except:
        trend[str(previous_chrm)+"-"+str(previous_binstart)]=QTLtrend
#        print(QTLtrend)

#%%
with open(str(out_dir)+str(name2)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+".csv","r") as file, open(str(out_dir)+str(name2)+"_QTL_regions_with_support_interval_LOD_"+str(args.lod)+"_with_trend.csv","w") as out:
    out.write("Chrm\tBin Start (bp)\tBin Stop (bp)\tPeak LOD\tPeak (bp)\tLOD threshold\tLOD support interval\tQTL trend\n")
    next(file)
    for line in file:
        ok = False
        line=line.split("\n")[0]
        line=line.split("\t")
#        print("\n",line)
        chrm = line[0]
        binstop = line[2]
#        print(binstop)
        for k,v in trend.items():
            k_chr = k.split("-")[0]
            k_binstop = k.split("-")[1]
            if k_binstop == binstop and k_chr == chrm:
                ok = True
#                print("  --YEEES",k,v)
                out.write("\t".join(line)+"\t"+"-".join(v)+"\n")
        if ok == False:
            print("  --NOT FOUND")
            break
#%%
