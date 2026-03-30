# Author Lvmin
# coding=utf-8
# @Time    : 2022/5/7 8:42
# @File    : count_intro_number.py
# @contact: lvmingjie_good@163.com
import re, sys, os, glob
file_path = sys.argv[1]
files = glob.glob(file_path+"/*continuous.fill.txt")
for intro_file in files:
    # intro_file = "E:/projects/people/weibo/2_weibo_wheat/introgression/2022May/matrix_table_XJ166_10M_2022May111424.xls"
    fin = open(intro_file, 'r')
    out_file = re.sub(r'txt$', 'stat.txt', intro_file)
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
            region = regions[i].split("_")[0]
            if region not in dic_num[sample]:
                tag = 0
                tag1 = 0
                if i == 0:
                    if dic[sample][i] == "1" and dic[sample][i+1] == "0":
                        tag += 1
                else:
                    if dic[sample][i] == "1" and dic[sample][i-1] == "0":
                        tag += 1
                if dic[sample][i] == "1":
                    tag1 += 1
            else:
                if i == 0:
                    if dic[sample][i] == "1" and dic[sample][i+1] == "0":
                        tag += 1
                else:
                    if dic[sample][i] == "1" and dic[sample][i-1] == "0":
                        tag += 1
                if dic[sample][i] == "1":
                    tag1 += 1
            dic_num[sample][region] = [tag, tag1]
            # print(tag)
    fot.write("Chrom")
    for i in range(len(samples)):
        # fot.write("\t"+samples[i]+"\t"+samples[i]+"_all")
        fot.write("\t" + samples[i])
    fot.write("\n")

    for region in dic_num[samples[0]]:
        fot.write(region)
        for sample in dic_num:
            # fot.write("\t"+str(dic_num[sample][region][0])+"\t"+str(dic_num[sample][region][1]))
            fot.write("\t" + str(dic_num[sample][region][0]))
        fot.write("\n")



