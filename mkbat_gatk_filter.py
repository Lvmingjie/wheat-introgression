# Author Lvmin
# coding=utf-8
# @Time    : 2023/11/17 15:42
# @File    : mkbat_gatk_filter.py
# @contact: lvmingjie_good@163.com
import glob, os, re
import argparse
parser = argparse.ArgumentParser(description='make bats for bwa and gatk.')
parser.add_argument('-r','--ref', dest='ref', metavar='', help='reference genome fasta file', required=True)  #默认4G内存
parser.add_argument('-d','--dir', dest='dir', type = str, metavar='', help='working root directory, the directory should at least have "1_raw or 1_raw/fastp" dir with *.R1.fastq.gz files', required=True)
# parser.add_argument('-N','--Node', dest='Node', metavar='1', help='number of cpu', required=False)
parser.add_argument('-m','--mem', dest='mem', metavar='50000', default="50G", help='memory xx Mb', required=False)  #默认10G内存
args = parser.parse_args()
CPU = '4'
GATK = '/public200T/apps/gatk-4.2.3.0/gatk'
BWA = '/public200T/apps/bwa/bwa'
SAMTOOLS = '/usr/local/bin/samtools'
REFER = args.ref
FAI = REFER + '.fai'
work_dir = args.dir
runbat = work_dir+"/run_bwa_gatk.bat"
fot = open(runbat, 'w')
chrs = []
fin_fai = open(FAI, 'r')
dic_seg = {}
win = 5000000

for line in fin_fai:
    chrom = line.strip().split("\t")[0]
    if re.match(r'^C|^c|^\d|^A|^D|^arahy.Tifrunner.gnm2.Arahy', chrom):
        chrs.append(chrom)
        chr_length = int(line.strip().split("\t")[1])
        for i in range(chr_length // win):
            if i == 0:
                start = 1
                end = (i+1)*win
                key = "%s:%s-%s" %(chrom, start, end)
                dic_seg[key] = [chrom, start, end]
            else:
                start = i*win
                end = (i+1)*win
                key = "%s:%s-%s" %(chrom, start, end)
                dic_seg[key] = [chrom, start, end]
        start = (chr_length // win) * win
        end = chr_length
        key = "%s:%s-%s" % (chrom, start, end)
        dic_seg[key] = [chrom, start, end]

VCF =  work_dir + '/4_vcf'
if not os.path.exists(VCF):
    os.mkdir(VCF)
if not os.path.exists(VCF+"/tmp"):
    os.mkdir(VCF+"/tmp")

VCF2 =  work_dir + '/5_vcf_filter'
if not os.path.exists(VCF2):
    os.mkdir(VCF2)

for window in dic_seg:
    chrom = window.split(":")[0]
    start = window.split(":")[1].split("-")[0]
    end = window.split(":")[1].split("-")[1]
    bat_snp = "%s/%s.filter_snp.bat" %(VCF2, window)
    bat_indel = "%s/%s.filter_indel.bat" %(VCF2, window)
    fot = open(bat_snp, 'w')
    fot1 = open(bat_indel, 'w')
    ###for snp
    sbatch_prefix_1 = "#!/bin/bash\n" \
                      "#SBATCH --ntasks-per-node=2\n" \
                      "#SBATCH --time=7-00:00:00\n" \
                      "#SBATCH --mem=20G\n" \
                      "#SBATCH --job-name=snp_%s\n" \
                      "#SBATCH --error=%s/snp_%s.err\n" \
                      "#SBATCH --output=%s/snp_%s.out\n\n\n" % ( window, VCF2, window, VCF2, window)
    cmd1 = "gatk SelectVariants -select-type SNP -V %s/Reg_%s_%s_%s.vcf.gz -O %s/Reg_%s_%s_%s.snp.vcf.gz \n" % (VCF, chrom, start, end, VCF2, chrom, start, end)
    cmd2 = 'gatk VariantFiltration -V %s/Reg_%s_%s_%s.snp.vcf.gz ' \
            '--filter-expression "QD < 2.0" --filter-name "QD2" ' \
            '--filter-expression "QUAL < 30.0" --filter-name "QUAL30" ' \
            '--filter-expression "SOR > 3.0" --filter-name "SOR3" ' \
            '--filter-expression "MQ < 40.0" --filter-name "MQ40.0" ' \
            '--filter-expression "FS > 60.0" --filter-name "FS60.0" ' \
            '--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" ' \
            '--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8.0" ' \
            '--cluster-size 3 --cluster-window-size 10 ' \
            '-O %s/Reg_%s_%s_%s.snp.filter.vcf.gz \n' % (VCF2, chrom, start, end, VCF2, chrom, start, end)
    cmd3 = 'bcftools +setGT %s/Reg_%s_%s_%s.snp.filter.vcf.gz -- -t q -n . -i "FMT/DP=0" |bcftools view -f PASS -Ov |' \
           'bcftools annotate -x  ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ' \
           ' -Oz -o %s/Reg_%s_%s_%s.snp.final.vcf.gz \n' %(VCF2, chrom, start, end, VCF2, chrom, start, end)
    # cmd3.1 = " bcftools annotate --set-id +'%CHROM\_%POS' -O z -o ${outPath}${File}_snp_filter"
    cmd4 = 'tabix %s/Reg_%s_%s_%s.snp.final.vcf.gz\n' %(VCF2, chrom, start, end)
    fot.write(sbatch_prefix_1+cmd1+cmd2+cmd3+cmd4)
    sbatch_prefix_2 = "#!/bin/bash\n" \
                      "#SBATCH --ntasks-per-node=2\n" \
                      "#SBATCH --time=7-00:00:00\n" \
                      "#SBATCH --mem=20G\n" \
                      "#SBATCH --job-name=indel_%s\n" \
                      "#SBATCH --error=%s/indel_%s.err\n" \
                      "#SBATCH --output=%s/indel_%s.out\n\n\n" % ( window, VCF2, window, VCF2, window)
    
    cmd5 = "gatk SelectVariants -select-type INDEL -V %s/Reg_%s_%s_%s.vcf.gz -O %s/Reg_%s_%s_%s.indel.vcf.gz \n" % (VCF, chrom, start, end, VCF2, chrom, start, end)
    cmd6 = 'gatk VariantFiltration -V %s/Reg_%s_%s_%s.indel.vcf.gz ' \
            '--filter-expression "QD < 2.0" --filter-name "QD2" ' \
            '--filter-expression "QUAL < 30.0" --filter-name "QUAL30" ' \
            '--filter-expression "FS > 200.0" --filter-name "FS200.0" ' \
            '--filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20.0" ' \
            '-O %s/Reg_%s_%s_%s.indel.filter.vcf.gz \n' % (VCF2, chrom, start, end, VCF2, chrom, start, end)
    cmd7 = 'bcftools +setGT %s/Reg_%s_%s_%s.indel.filter.vcf.gz -- -t q -n . -i "FMT/DP=0" |bcftools view -f PASS -Ov |' \
           'bcftools annotate -x  ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ' \
           ' -Oz -o %s/Reg_%s_%s_%s.indel.final.vcf.gz \n' %(VCF2, chrom, start, end, VCF2, chrom, start, end)
    cmd8 = 'tabix %s/Reg_%s_%s_%s.indel.final.vcf.gz\n' %(VCF2, chrom, start, end)
    fot1.write(sbatch_prefix_2+cmd5+cmd6+cmd7+cmd8)
