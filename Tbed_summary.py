#!/usr/bin/env python3

import os,sys,re,pysam

def xopen(bedfile,genes):
    genes = [i.upper() for i in genes]
    genes = set(genes)
    #bed_dic = {}
    #i = 0
    with open(bedfile) as handle:
        for line in handle:
            if line:
                items = line.strip().split(None)

                if len(items) != 4: continue
                chrom = str(items[0])
                start = int(items[1])
                end = int(items[2])
                gene_ex = items[3]
                #gene_bed = [i.split(".")[0].upper() for i in gene_ex.split(";")]
                gene_bed = [i.upper() for i in gene_ex.split(";")]
                intersect=list(set(gene_bed) & genes)
                if len(intersect) == 0 : continue
                #i += 1
                #bed_dic.update({i:[chrom,start,end,gene_ex]})
                yield chrom, start, end, gene_ex
        #return bed_dic

def tackle_bam(bamfile,chrom,start,end,gene_ex,all_region_reads,all_pos_cover,pos_cont_sum):
    samfile = pysam.AlignmentFile(bamfile,'rb')

    region_reads = samfile.count(chrom,start,end)    #result same with samfile.fetch
    all_region_reads += region_reads    

    i = 0
    #for read in samfile.fetch(chrom,start,end):
    #    i += 1
    #print(chrom,start,end,i,"samfile.fetch")   #result same with samfile.count

    pos_cont = []
    for pc in samfile.pileup(chrom,start,end,max_depth=900000,truncate=True):
        pos_cont.append(pc.n)
        #if pc.n < 501:
        #   print ("pos:", pc.pos, pc.n, file=sys.stderr)
        #print("POS",chrom,start,end)
        #for pcread in pc.pileups:
        #    print(pcread)
    samfile.close()

    pos_full = ':'.join([chrom,str(start + 1),str(end)])
    all_pos_cover.update({pos_full:pos_cont})
    pos_cont_sum += sum(pos_cont)

    return pos_cont_sum,all_region_reads,all_pos_cover

def onebed(bedfile,bamfile,genes):
    bed_length = 0;all_base = 0;all_region_reads = 0; pos_cont_sum = 0
    all_pos_cover = {}
    for chrom, start, end, gene_ex  in  xopen(bedfile,genes):
        bed_length += end - start + 1
        pos_cont_sum,all_region_reads,all_pos_cover = tackle_bam(bamfile,chrom,start,end,gene_ex,all_region_reads,all_pos_cover,pos_cont_sum)
    #print("pos_cont_sum",pos_cont_sum)
    #print("all_region_reads",all_region_reads)
    #print("all_pos_cover",all_pos_cover)
    #print("bed_length",bed_length)
    return pos_cont_sum,all_region_reads,all_pos_cover,bed_length   

def get_all_reads(bedfile,bamfile,genes):
    genes = '|'.join(genes)
    ID = os.path.basename(bamfile).split(".")[0]
    bed_cmd = 'grep -iP "{genes}" {bedfile} > {ID}.bed'.format(genes=genes,bedfile=bedfile,ID=ID)
    os.system(bed_cmd)
    all_cmd = 'samtools view -L {ID}.bed {bamfile} > {ID}_all_tmd_reads.sam'.format(bamfile=bamfile,ID=ID,)
    mapped_cmd = 'samtools view -F 4 -L {ID}.bed {bamfile} > {ID}_mapped_tmd_reads.sam'.format(bamfile=bamfile,ID=ID)
    os.system(mapped_cmd)
    os.system(all_cmd)
    mapped_file = 'wc -l {ID}_mapped_tmd_reads.sam'.format(ID=ID)
    all_file = 'wc -l {ID}_all_tmd_reads.sam'.format(ID=ID)

    all_reads = int(list(os.popen(all_file))[0].split(None)[0].strip())   #r1
    mapped_reads = int(list(os.popen(mapped_file))[0].split(None)[0].strip()) #r2
    pt_reads = round(mapped_reads/all_reads * 100,2) #r3
    #print("----all_reads",all_reads)
    all_bases = all_reads * 150 #r4
    #print("----all_bases",all_bases)
    soft_bases = 'grep -oP "\dS" {ID}_all_tmd_reads.sam'.format(ID=ID)
    #print("-----Soft",soft_bases)
    soft_bases = [int(i.strip().replace("S","")) for i in list(os.popen(soft_bases))]
    all_target_read = all_bases - sum(soft_bases) #r5
    #print("----all_target_read",all_target_read)
    pt_base = round(all_target_read/all_bases * 100,2) #r6

    var_bases = 'grep -oP "NM:i:[0-9]*" {ID}_all_tmd_reads.sam'.format(ID=ID)
    var_bases = [int(i.strip().split(":")[2]) for i in list(os.popen(var_bases))]
    all_target_bases = all_bases - sum(var_bases)
    #print("----all_target_bases",all_target_bases)
    pt_concor = round(all_target_bases/all_bases * 100,2) #7


    os.remove('{ID}_all_tmd_reads.sam'.format(ID=ID))
    os.remove('{ID}_mapped_tmd_reads.sam'.format(ID=ID))
    os.remove('{ID}.bed'.format(ID=ID))
  
    return all_reads,mapped_reads,pt_reads,all_bases,all_target_read,pt_base,pt_concor 

    
 
def main(bedfile_real,bamfile,genes):
    all_reads,mapped_reads,pt_reads,all_bases,all_target_read,pt_base,pt_concor = get_all_reads(bedfile_real,bamfile,genes)
    #print("all_reads",all_reads)
    #print("mapped_reads",mapped_reads)
    #print("pt_reads",pt_reads)
    #print("all_bases",all_bases)
    #print("all_target_read",all_target_read)
    #print("pt_base",pt_base)
    #print("pt_concor",pt_concor)

    real_pos_cont,real_reads_cont,real_pos_cover,real_bed_length = onebed(bedfile_real,bamfile,genes)
    #print("target_length",real_bed_length)
    #print("average_depth",round(real_pos_cont/real_bed_length,0))
    bases = 0;X1 = 0;X20 = 0; X100 = 0; X500 = 0
    for key in real_pos_cover:
        for base_dp in real_pos_cover[key]:
            bases += 1
            if base_dp >= 1:
                X1 += 1
            if base_dp >= 20:
                X20 += 1
            if base_dp >= 100:
                X100 += 1
            if base_dp >= 500:
                X500 += 1
    print("all_reads",all_reads)
    print("pt_reads",str(pt_reads)+"%")
    print("all_bases",all_bases)
    print("all_target_read",all_target_read)
    print("pt_base",str(pt_base)+"%")
    print("target_length",real_bed_length)
    print("average_depth",round(real_pos_cont/real_bed_length,0))
    print("pt_concor",str(pt_concor)+"%")
    #print("bases",bases,real_pos_cont)
    base = real_bed_length
    print(">=1X",str(round(X1/bases * 100,2))+"%")
    print(">=20X",str(round(X20/bases * 100,2))+"%")
    print(">=100X",str(round(X100/bases * 100,2))+"%")
    print(">=500X",str(round(X500/bases * 100,2))+"%")

if __name__ == "__main__":
    try:
        bedfile_real = sys.argv[1]
        bamfile = sys.argv[2]
        genes = sys.argv[3:]
        main(bedfile_real,bamfile,genes)
    except:
        print("example: ","Tbed_summary.py oncosmart2.b37.bed S170000999-PRL.realn.bam ALK BRAF EGFR ERBB2 KIT KRAS MET PIK3CA RET ROS1")
