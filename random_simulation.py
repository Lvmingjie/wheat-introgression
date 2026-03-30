# Author Lvmin
# coding=utf-8
# @Time    : 2022/7/14 14:32
# @File    : random_simulation.py
# @contact: lvmingjie_good@163.com
import re, glob, sys, random, linecache, os
file_path = sys.argv[1]
files = glob.glob(file_path+"/*_1M*.continuous.fill.txt")
for file in files:
    print(file)
    key = re.match(r'.+_(.+\d+)_1M', os.path.basename(file)).group(1)[:2]
    print(key)
    fot = open(file_path+"/simulation/simulation."+key+".txt", 'w')
    fin = open(file, 'r')
    ti = fin.readline().strip("\n")
    samples = ti.split("\t")[4:]

    dic = {}
    regions = []
    for line in fin:
        line = line.strip("\n")
        _line = line.split("\t")
        for i in range(len(_line[4:])):
            if samples[i] not in dic:
                dic[samples[i]] = {}
                dic[samples[i]][_line[0]] = _line[4 + i]
            else:
                dic[samples[i]][_line[0]] = _line[4 + i]
        regions.append(_line[0])
    # print(dic["XJ_4500_01"])
    # print(samples)
    rdm_samples = random.sample(samples, 10)
    # print(rdm_samples)
    i = 1
    select_samples = []
    dic_num = {}
    while i <= len(samples):
        dic_num[i] = {}
        remain_samples = [x for x in samples if x not in select_samples]
        get_one = random.sample(remain_samples, 1)
        select_samples = select_samples + get_one
        # print(select_samples)
        for sample in select_samples:
            # dic_num[i][sample] = {}
            for region in dic[sample]:
                if dic[sample][region] == "0":
                    continue
                else:
                    dic_num[i][region] = "1"
        i += 1

    for i in range(1, len(samples)+1):
        pct = len(dic_num[i])/len(regions)
        fot.write("%s\t%.4f\n" % (i, pct))
