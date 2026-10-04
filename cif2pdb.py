#!/usr/bin/env python3

import sys 
from math import *
import datetime
import gzip
import re

def preprocess(lines,filename):
    for line in lines:
        if line[0] == "#":
            continue

        if line[0:5] == 'data_':
            break
        else:
            print("{} is not mmCIF data".format(filename))
            exit()

    newlines = []
    dummy = ''

    total = len(lines)
    i = 0
    while i < total:
        if len(lines[i]) == 0:
            i += 1
            continue
        
        if lines[i][0] == ";":
            dummy += '"'
            dummy += (lines[i][1:].replace("'","").replace('"',''))
            i += 1
            while i < total and len(lines[i]) == 0:
                i += 1
            while i < total and lines[i][0] != ";":
                dummy += (lines[i].replace("'","").replace('"',''))
                i += 1
                while i < total and len(lines[i]) == 0:
                    i += 1
            dummy += '"'
            newlines.append(dummy)
            dummy = ''
        else:
            newlines.append(lines[i])
        i += 1

    return newlines
            
def writeHEADER(pdb,filename,f):
    f.write("HEADER    CONVERTED FROM CIF TO PDB FORMAT        XX-XXX-XX   {}\n".format(pdb['_entry.id']))

    if '_entity.pdbx_description' in pdb:
        if type(pdb['_entity.pdbx_description']) is list:
            dscrpt = ' / '.join(pdb['_entity.pdbx_description']).replace(" / water","").replace("'","")
        elif type(pdb['_entity.pdbx_description']) is str:
            dscrpt = pdb['_entity.pdbx_description'].replace("'","")
        else:
            dscrpt = pdb['_entity.pdbx_description']
        f.write("REMARK    ")
        for i in range(len(dscrpt)):
            if (i+1) % 70 == 0:
                f.write("\nREMARK    ")
            f.write(dscrpt[i])
        f.write("\n")
    else:
        f.write("REMARK    No _entity.pdbx_description exists.\n")

    f.write("REMARK    {}\n".format(datetime.datetime.now()))
    f.write("REMARK    CONVERTED FROM {}\n" .format(filename))
    return

def writeATOM(pdb,f):
    ln = 1
    total = len(pdb['_atom_site.group_PDB'])
    cnt = 0
    old_asym_id = pdb['_atom_site.label_asym_id'][0]
    for i in range(total):            
        if pdb['_atom_site.group_PDB'][i] == "ATOM":
            if '_atom_site.pdbx_PDB_model_num' in pdb:
                if cnt == 0:
                    model = pdb['_atom_site.pdbx_PDB_model_num'][i]
                    cnt += 1
                elif pdb['_atom_site.pdbx_PDB_model_num'][i] != model:
                    continue
            
            if old_asym_id != pdb['_atom_site.label_asym_id'][i]:
                f.write("TER  \n")
                old_asym_id = pdb['_atom_site.label_asym_id'][i]

            if log10(int(pdb['_atom_site.id'][i])) >= 7.0:
                f.write("ATOM 999999")
                print("ATOM: atom number, compromised",file=sys.stderr)
            else:
                f.write("ATOM{0:7d}".format(int(pdb['_atom_site.id'][i])))
                
            l = len(pdb['_atom_site.label_atom_id'][i])
            if l == 1:
                f.write("  {0}   {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                                pdb['_atom_site.label_comp_id'][i]))
            elif l == 2:
                f.write("  {0}  {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                               pdb['_atom_site.label_comp_id'][i]))
            elif l == 3:
                f.write("  {0} {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                              pdb['_atom_site.label_comp_id'][i]))
            elif l == 4:
                f.write(" {0} {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                             pdb['_atom_site.label_comp_id'][i]))
            else:
                f.write("{0} {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                            pdb['_atom_site.label_comp_id'][i]))
                
            if len(pdb['_atom_site.label_asym_id'][i]) == 1:
                if pdb['_atom_site.label_seq_id'][i].isdecimal() == False:
                    f.write(" {0:1s}{1:4d}".format(pdb['_atom_site.label_asym_id'][i], \
                                            int('1')))
                else:
                    f.write(" {0:1s}{1:4d}".format(pdb['_atom_site.label_asym_id'][i], \
                                            int(pdb['_atom_site.label_seq_id'][i])))
            elif len(pdb['_atom_site.label_asym_id'][i]) == 2:
                if int(pdb['_atom_site.label_seq_id'][i]) >= 1000:
                    num = int(pdb['_atom_site.label_seq_id'][i][-3:])
                    print("ATOM: residue number, compromised",file=sys.stderr)
                else:
                    num = int(pdb['_atom_site.label_seq_id'][i])                    
                f.write(" {0:2s}{1:3d}".format(pdb['_atom_site.label_asym_id'][i],num))
            else:
                if int(pdb['_atom_site.label_seq_id'][i]) >= 1000:
                    num = int(pdb['_atom_site.label_seq_id'][i][-3:])
                    print("ATOM: residue number, compromised",file=sys.stderr)
                else:
                    num = int(pdb['_atom_site.label_seq_id'][i])
                f.write(" {0:2s}{1:3d}".format(pdb['_atom_site.label_asym_id'][i][0:2],num))
                print("ATOM: chain id, compromised",file=sys.stderr)
                
            if pdb['_atom_site.label_seq_id'][i].isdecimal() == False:
                ln = 1 
            else:
                ln = int(pdb['_atom_site.label_seq_id'][i])

            if pdb['_atom_site.label_alt_id'][i] == '.':
                f.write(" ")
            else:
                f.write("{0:1s}".format(pdb['_atom_site.label_alt_id'][i]))

            f.write("   {x:8.3f}{y:8.3f}{z:8.3f}{o:6.2f}{b:6.2f}{s:>12s}\n".format(\
                      x=float(pdb['_atom_site.Cartn_x'][i]),\
                      y=float(pdb['_atom_site.Cartn_y'][i]),\
                      z=float(pdb['_atom_site.Cartn_z'][i]),\
                      o=float(pdb['_atom_site.occupancy'][i]),\
                      b=float(pdb['_atom_site.B_iso_or_equiv'][i]),\
                      s=pdb['_atom_site.type_symbol'][i]))
           
    f.write("TER  \n")        
    return ln

def writeHETATM(pdb,ln,f):
    total = len(pdb['_atom_site.group_PDB'])
    for i in range(total):
        if pdb['_atom_site.group_PDB'][i] == "HETATM":
            if '_atom_site.pdbx_PDB_model_num' in pdb:
                if pdb['_atom_site.pdbx_PDB_model_num'][i] != '1':
                    continue

            if log10(int(pdb['_atom_site.id'][i])) >= 5.0:
                f.write("HETATM99999")
                print("HETATM: atom number, compromised",file=sys.stderr)
            else:
                f.write("HETATM{0:5d}".format(int(pdb['_atom_site.id'][i])))
                
            l = len(pdb['_atom_site.label_atom_id'][i])
            m = min(3,len(pdb['_atom_site.label_comp_id'][i]))

            if l == 1:
                f.write("  {0}   {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                                pdb['_atom_site.label_comp_id'][i][0:m]))
            elif l == 2:
                f.write("  {0}  {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                               pdb['_atom_site.label_comp_id'][i][0:m]))
            elif l == 3:
                f.write("  {0} {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                              pdb['_atom_site.label_comp_id'][i][0:m]))
            elif l == 4:
                f.write(" {0} {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                             pdb['_atom_site.label_comp_id'][i][0:m]))
            else:
                f.write("{0} {1:>3s}".format(pdb['_atom_site.label_atom_id'][i], \
                                            pdb['_atom_site.label_comp_id'][i][0:m]))

            if pdb['_atom_site.label_seq_id'][i] == '.':
                if pdb['_atom_site.label_asym_id'][i-1] != pdb['_atom_site.label_asym_id'][i]:
                    ln += 1
                    if ln >= 1000:
                        ln = 1 
                        print("HETATM: residue number, compromised",file=sys.stderr)
                if len(pdb['_atom_site.label_asym_id'][i]) == 1:
                    f.write(" {0:1s}{1:4d}".format(pdb['_atom_site.label_asym_id'][i],ln))
                elif len(pdb['_atom_site.label_asym_id'][i]) == 2:
                    f.write(" {0:2s}{1:3d}".format(pdb['_atom_site.label_asym_id'][i],ln))
                else:
                    f.write(" {0:2s}{1:3d}".format(pdb['_atom_site.label_asym_id'][i][0:2],ln))
                    print("HETATM: chain id, compromised",file=sys.stderr)
            else:
                if len(pdb['_atom_site.label_asym_id'][i]) == 1:
                    f.write(" {0:1s}{1:4d}".format(pdb['_atom_site.label_asym_id'][i],\
                                               int(pdb['_atom_site.label_seq_id'][i])))
                elif len(pdb['_atom_site.label_asym_id'][i]) == 2:
                    if int(pdb['_atom_site.label_seq_id'][i]) >= 1000:
                        num = int(pdb['_atom_site.label_seq_id'][i][-3:])
                        print("HETATM: residue number, compromised (2)",file=sys.stderr)
                    else:
                        num = int(pdb['_atom_site.label_seq_id'][i])
                    f.write(" {0:2s}{1:3d}".format(pdb['_atom_site.label_asym_id'][i],num)) 
                else:
                    f.write(" {0:2s}{1:3d}".format(pdb['_atom_site.label_asym_id'][i][0:2],num)) 
                    print("HETATM: chain id, compromised (2)",file=sys.stderr)
                    
            if pdb['_atom_site.label_alt_id'][i] == '.':
                f.write(" ")
            else:
                f.write("{0:1s}".format(pdb['_atom_site.label_alt_id'][i]))

            f.write("   {x:8.3f}{y:8.3f}{z:8.3f}{o:6.2f}{b:6.2f}{s:>12s}\n".format(\
                       x=float(pdb['_atom_site.Cartn_x'][i]),\
                       y=float(pdb['_atom_site.Cartn_y'][i]),\
                       z=float(pdb['_atom_site.Cartn_z'][i]),\
                       o=float(pdb['_atom_site.occupancy'][i]),\
                       b=float(pdb['_atom_site.B_iso_or_equiv'][i]),\
                       s=pdb['_atom_site.type_symbol'][i]))
    return

def writeEND(pdb,f):
    f.write("END  \n")
    return

def readCIF(lines,pdb):
    l = len(lines)
    i = 0
    while i < l:
        if lines[i][0] == "#":
            i += 1
            continue

        if lines[i][0:5] == "data_":
            print(lines[i],file=sys.stderr)
            i += 1
            continue
        
        if lines[i][0] == "_":
            el = lines[i].split()
            pdb[el[0]] = ' '.join(el[1:])
            if len(pdb[el[0]]) == 0:
                i += 1
                while lines[i][0] != '_' and lines[i][0] != '#' and "loop_" not in lines[i]:
                    el2 = lines[i].split()
                    pdb[el[0]] += ' '.join(el2)
                    i += 1
                i -= 1
        elif lines[i][0:5] == "loop_":
            i += 1
            keys = list()
            while lines[i][0] =='_':
                keys.append(lines[i].split()[0])
                pdb[keys[-1]] = []
                i += 1
            while i < l and lines[i][0] != '_' and lines[i][0] != '#' and "loop_" not in lines[i]:
                el = re.findall(r"'[^']*'|\S+", lines[i])
                if lines[i].count('"') % 2 != 0:
                     lines[i] += '"'
                if lines[i].count("'") % 2 != 0:
                     lines[i] += "'"
                el = re.findall(r"'[^']*'|\S+", lines[i])
                for j, key in enumerate(keys):
                    if i+1 < l and len(el) <= j:
                        i += 1
                        while i+1 < l and lines[i][0] == '#':
                            i += 1
                        if i < l and (lines[i][0] == '_' or "loop_" in lines[i]):
                            i -= 1
                            break
                        if i < l:
                            el.extend(re.findall(r"'[^']*'|\S+", lines[i]))
                            if lines[i].count('"') % 2 != 0:
                                lines[i] += '"'
                            if lines[i].count("'") % 2 != 0:
                                lines[i] += "'"
                            el.extend(re.findall(r"'[^']*'|\S+", lines[i]))
                    if j < len(el):
                        pdb[key].append(el[j].replace('"',''))
                i += 1
            i -= 1
        else:
            print("Irregular line: ",lines[i],file=sys.stderr)
            
        i += 1
        
    return
#
#
if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("%s   cif    output" % sys.argv[0],file=sys.stderr)
        exit ()
        
    pdb = {}
    suffix = ".gz"
    if sys.argv[1].endswith(suffix):
        with gzip.open(sys.argv[1],'rt') as f:
            alllines = f.read().split("\n")
            alllines.pop()    
    else:
        with open(sys.argv[1],'r') as f:
            alllines = f.read().split("\n")
            alllines.pop()

    alllines = preprocess(alllines,sys.argv[1])
    readCIF(alllines,pdb)

    with open(sys.argv[2],"x") as f:
        writeHEADER(pdb,sys.argv[1],f)
        ln = writeATOM(pdb,f)
        writeHETATM(pdb,ln,f)
        writeEND(pdb,f)
#EOF
