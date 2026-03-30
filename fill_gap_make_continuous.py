# Author Lvmin
# coding=utf-8
# @Time    : 2022/5/10 14:36
# @File    : make_continuous_0.py
# @contact: lvmingjie_good@163.com
import re, sys, os, glob
file_path = sys.argv[1]
files = glob.glob(file_path+"/*_1M*continuous.txt")
for intro_file in files:
    print(intro_file)
    fin = open(intro_file, 'r')
    out_file = re.sub(r'continuous.txt$', 'continuous.fill.txt', intro_file)
    fot = open(out_file, 'w')
    ti = fin.readline().strip("\n")
    fot.write(ti+"\n")
    _ti = ti.split("\t")
    samples = _ti[2:]
    dic = {}
    ld = 5   #LD
    tag = 0
    regions = {}
    color = {}
    for line in fin:
        line = line.strip("\n")
        _line = line.split("\t")
        for i in range(len(_line[2:])):
            value = _line[2+i]
            if samples[i] not in dic:
                dic[samples[i]] = {}
                dic[samples[i]][tag] = _line[2+i]
            else:
                dic[samples[i]][tag] = _line[2+i]
        regions[tag] = _line[0]
        tag += 1
        color[_line[0]] = _line[1]

    dic_gaps = {}
    for sample in samples:
        # print(sample)
        for i in range(tag):
            if sample not in dic_gaps:
                dic_gaps[sample] = []
                if dic[sample][i] == "0":
                    gaps = [i]
                else:
                    gaps = []
            else:
                if dic[sample][i] == "0":
                    gaps.append(i)
                    # print(gaps)
                else:
                    # print(dic_gaps[sample])
                    if len(gaps) >= 1:
                        dic_gaps[sample].append(gaps)
                    gaps = []
        dic_gaps[sample].append(gaps)
    dic1 = {}
    for sample in samples:
        dic1[sample] = {}
        for i in range(tag):
            dic1[sample][i] = dic[sample][i]

    for sample in samples:
        for i in range(1,len(dic_gaps[sample])-1):  # 
            if len(dic_gaps[sample][i]) <= ld:      #
                for j in dic_gaps[sample][i]:
                    dic1[sample][j] = "1"


    for i in range(len(dic[samples[0]])):
        fot.write(regions[i]+"\t"+color[regions[i]])
        for j in range(len(samples)):
            fot.write("\t"+dic1[samples[j]][i])
        fot.write("\n")
