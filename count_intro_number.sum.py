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
    out_file = re.sub(r'txt$', 'stat.sum.txt', intro_file)
    fot = open(out_file, 'w')
    ti = fin.readline().strip("\n")
    _ti = ti.split("\t")
    samples = _ti[4:]
    dic = {}
    regions = []
    dic_all = {}
    j = 0
    print(intro_file)
    for line in fin:
        line = line.strip("\n")
        _line = line.split("\t")
        dic_all[j] = "0"
        for i in range(len(_line[4:])):
            if samples[i] not in dic:
                dic[samples[i]] = [_line[4+i]]
            else:
                dic[samples[i]].append(_line[4+i])
            if _line[4+i] == "1":
                dic_all[j] = "1"  #unique
        j += 1
        # print(j)
        regions.append(_line[0])
    # print(dic_all[13])
    dic_num_seg_all = {}  #total segments
    # dic_num_seg_uniq = {}  #total segments
    dic_len_all = {}   #total introgression number
    dic_len_uniq = {}  #total unique introgression
    dic_chr_len = {}  #total chromosome length
    for i in range(len(dic_all)):
        region = regions[i].split("_")[0]
        if region not in dic_num_seg_all:
            chr_len = 1
            tag_all = 0   #segment number
            tag_all2 = 0  #intro legnth (1Mb)
            for sample in dic:
                if i == 0:
                    if dic[sample][i] == "1" and dic[sample][i+1] == "0":
                        tag_all += 1
                else:
                    if dic[sample][i] == "1" and dic[sample][i-1] == "0":
                        tag_all += 1
                if dic[sample][i] == "1":
                    tag_all2 += 1
        else:
            chr_len += 1
            for sample in dic:
                if i == 0:
                    if dic[sample][i] == "1" and dic[sample][i+1] == "0":
                        tag_all += 1
                else:
                    if dic[sample][i] == "1" and dic[sample][i-1] == "0":
                        tag_all += 1
                if dic[sample][i] == "1":
                    tag_all2 += 1

        # print(i, tag_all)
        dic_num_seg_all[region] = tag_all  #introgression number
        dic_len_all[region] = tag_all2
        dic_chr_len[region] = chr_len   #all length

    for i in range(len(dic_all)):
        region = regions[i].split("_")[0]
        if region not in dic_len_uniq:
            tag_all1 = 0
            if dic_all[i] == "1":
                tag_all1 += 1
        else:
            if dic_all[i] == "1":
                tag_all1 += 1
        dic_len_uniq[region] = tag_all1  #unique length

    fot.write("Chrom\tTotal number of segments in population\tAverage number of segments per sample\tTotal length of segments\tAverage length of segments\tTotal unique length of introgression\tPopulation coverage rate\n")


    for region in dic_num_seg_all:
        # print(region)
        fot.write(region+"\t"+str(dic_num_seg_all[region])+"\t"+str(dic_num_seg_all[region]/len(samples))+"\t"+str(dic_len_all[region])+"\t"+str(dic_len_all[region]/dic_num_seg_all[region])+"\t"+str(dic_len_uniq[region])+"\t"+str(dic_len_uniq[region]/dic_chr_len[region])+"\n")



