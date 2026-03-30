# Author Lvmin
# coding=utf-8
# @Time    : 2022/5/12 14:22
# @File    : count_intro_length.py
# @contact: lvmingjie_good@163.com
import re, sys, os, glob
file_path = sys.argv[1]
files = glob.glob(file_path+"/*continuous.txt")
for intro_file in files:
    fin = open(intro_file, 'r')
    out_file = re.sub(r'txt$', 'stat_len.txt', intro_file)
    fot = open(out_file, 'w')
    ti = fin.readline().strip("\n")
    _ti = ti.split("\t")
    samples = _ti[4:]
    dic = {}
    regions = []
    for line in fin:
        line = line.strip("\n")
        _line = line.split("\t")
        for i in range(len(_line[4:])):
            if samples[i] not in dic:
                dic[samples[i]] = [_line[4+i]]
            else:
                dic[samples[i]].append(_line[4+i])
        regions.append(_line[0])

    dic_num = {}
    for sample in dic:
        dic_num[sample] = {}
        for i in range(len(dic[sample])):
            if i == 0:
                if dic[sample][i] == "1" and dic[sample][i+1] == "0":
                    dic_num[sample][i] = 1
                elif dic[sample][i] == "1" and dic[sample][i+1] == "1":
                    tag = 1
                else:
                    continue
            else:
                if dic[sample][i] == "1":
                    if i == len(dic[sample])-1:
                        if dic[sample][i-1] == "0":
                            dic_num[sample][i] = 1
                        else:
                            tag += 1
                            dic_num[sample][i - tag + 1] = tag
                    else:
                        if dic[sample][i-1] == "0" and dic[sample][i+1] == "0":
                            dic_num[sample][i] = 1
                        elif dic[sample][i-1] == "0":
                            tag = 1
                        elif dic[sample][i+1] == "0":
                            tag += 1
                            dic_num[sample][i-tag+1] = tag
                        else:
                            tag += 1
                else:
                    continue


    for i in range(len(samples)):
        for j in dic_num[samples[i]]:
            fot.write(samples[i].split("_")[0]+"\t"+samples[i]+"\t"+regions[j].split("_")[0]+"\t"+str(dic_num[samples[i]][j])+"\n")
