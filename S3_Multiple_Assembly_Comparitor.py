from Bio import SeqIO
import os, sys

def ORFs_predictor(assembly):
    try:
        print 'Predicting ORFs of', assembly
        os.system('prodigal -a temp.orfs.faa -i '+str(assembly)+' -d temp.orfs.nt.faa -m -o temp.txt -p meta -q')
        os.system('cut -f1 -d \" \" temp.orfs.faa > '+str(assembly)+'.orfs.faa')
        os.system('cut -f1 -d \" \" temp.orfs.nt.faa > '+str(assembly)+'.orfs.fna')
    except:
        print 'ORFs prediction error! Please make sure prodigal is installed in your system'
    return str(assembly)+'.orfs.fna'

def ORFs_aligner(ORFs_assembly1, ORFs_assembly2):
    print 'Using BLAST to align', ORFs_assembly2, 'to', ORFs_assembly1
    try:
        blast_output=str(ORFs_assembly2)+'_vs_'+str(ORFs_assembly1)+'.txt'
        os.system('makeblastdb -in '+str(ORFs_assembly1)+' -dbtype nucl -hash_index -parse_seqids')
        os.system('blastn -query '+str(ORFs_assembly2)+' -db '+str(ORFs_assembly1)+' -evalue 1e-20 -max_target_seqs 10 -outfmt 6 -num_threads 42 -out '+str(blast_output))
    except:
        print 'Alignment error! Please check whether BLAST+ is installed in your system'
    return str(blast_output)

def ORFs_sequence_length_recorder(ORF):
    print 'Reading', ORF, 'ORFs length'
    seq_len={}
    for record in SeqIO.parse(ORF,'fasta'):
        seq_len[str(record.id)]=len(record.seq)
    print '---------------------' 
    return seq_len

def genome_contigs_recorder(binset, binset_record, binset_genome_size, coverage_matrix):
    binset_record[str(binset)]={}
    binset_genome_size[str(binset)]={}
    bins_coverage, bin_num_contigs, bin_total_length, bin_GC, bin_GC_ratio, binset_coverage, all_bins={}, {}, {}, {}, {}, {}, {}
    pwd=os.getcwd()

    print 'Parsing', binset

    n1=0
    for line in open(coverage_matrix, 'r'):
        n1+=1
        if n1 == 1:
            num=str(line).strip().count('drange')
        else:
            ids=str(line).strip().split('\t')[0]
            bins_coverage[str(ids)]={}
            i=1
            while i <= int(num):
                bins_coverage[str(ids)][i]=str(line).strip().split('\t')[3*i+1]
                i+=1

    for root, dirs, files in os.walk(binset):
        os.chdir(pwd+'/'+binset)
        for file in files:
            if '_genomes.' in file:
                hz=file.split('_genomes.')[1]
                if '.fasta' in hz or '.fa' in hz:
                    binset_coverage[file]={}
                    bin_total_length[file]=0
                    bin_GC[file]=0
                    n=0
                    i=1
                    while i <= int(num):
                        binset_coverage[file][i]=0
                        i+=1

                    for record in SeqIO.parse(file,'fasta'):        
                        n+=1
                        bin_num_contigs[file]=n
                        bin_total_length[file]+=len(record.seq)
                        bin_GC[file]+=int(str(record.seq).count('G'))
                        bin_GC[file]+=int(str(record.seq).count('C'))
                        if n == 1:
                            binset_genome_size[str(binset)][str(file)]=len(record.seq)
                        else:
                            binset_genome_size[str(binset)][str(file)]+=len(record.seq)

                        if str(record.id) not in binset_record[str(binset)].keys():
                            binset_record[str(binset)][str(record.id)]=str(file)
                        else:
                            binset_record[str(binset)][str(record.id)]+=','+str(file)
                    
                        i=1
                        while i <= int(num):
                            binset_coverage[file][i]+=float(bins_coverage[str(record.id)][i])
                            i+=1
                
                    if 'noclass' not in file and 'unbined' not in file:
                        all_bins[file]=1

    os.chdir(pwd)

    # print binset_coverage
    binset_coverage_avg={}
    f=open(str(binset)+'_bins_coverage.txt','w')
    i=1
    a='Bin'
    while i <= int(num):
        a+='\t'+'Coverage'+str(i)
        i+=1
    f.write(str(a)+'\n')

    for item in binset_coverage.keys():
        binset_coverage_avg[item]={}
        a=str(item)
        i=1
        while i <= int(num):
            avg_coverage=float(binset_coverage[item][i])/int(bin_num_contigs[item])
            binset_coverage_avg[item][i]=avg_coverage
            a+='\t'+str(avg_coverage)
            i+=1
        f.write(str(a)+'\n')
    f.close()

    f_t=open(binset+'_gc.txt', 'w')
    for item in bin_GC.keys():
        bin_GC_ratio[item]=round(100*float(bin_GC[item])/float(bin_total_length[item]),1)
        f_t.write(str(item)+'\t'+str(bin_GC_ratio[item])+'\n')
    f_t.close()
    print '---------------------' 
    
    return binset_record, binset_genome_size, binset_coverage_avg, bin_GC_ratio, all_bins

def coverage_GC_comparitor(binset_coverage_avg1, binset_coverage_avg2, binset_GC_ratio1, binset_GC_ratio2, iteration_num):
    print 'Comparing coverages'
    bins_score, bins_score_delta, bin_gc={}, {}, {}
    for item in binset_coverage_avg1.keys():
        num=len(binset_coverage_avg1[item])
        i=1
        bins_score[item], bins_score_delta[item]={}, {}
        while i <= num:
            a=str(binset_coverage_avg1[item][i])
            for item2 in binset_coverage_avg2.keys():
                b=str(binset_coverage_avg2[item2][i])
                ave_cov=(float(a)+float(b))/2
                delta_coverage=abs(float(b)-float(a))
                if float(ave_cov) != 0:
                    delta_coverage_vari=round(100*float(delta_coverage)/float(ave_cov), 2) # lower the delta_coverage is, closer the two bins
                else:
                    delta_coverage_vari=round(100*float(delta_coverage)/0.001, 2) # lower the delta_coverage is, closer the two bins
                    
                if ave_cov <= 1:
                    if item2 not in bins_score[item].keys():
                        bins_score[item][item2]=1
                        bins_score_delta[item][item2]=float(delta_coverage_vari) 
                    else:
                        bins_score[item][item2]+=1
                        bins_score_delta[item][item2]+=float(delta_coverage_vari)
                elif ave_cov > 1 and ave_cov <= 5:
                    if delta_coverage_vari <= 60: ### variable percentage
                        if item2 not in bins_score[item].keys():
                            bins_score[item][item2]=1
                            bins_score_delta[item][item2]=float(delta_coverage_vari) 
                        else:
                            bins_score[item][item2]+=1
                            bins_score_delta[item][item2]+=float(delta_coverage_vari)
                elif ave_cov > 5 and ave_cov <= 15:
                    if delta_coverage_vari <= 30: ### variable percentage
                        if item2 not in bins_score[item].keys():
                            bins_score[item][item2]=1
                            bins_score_delta[item][item2]=float(delta_coverage_vari) 
                        else:
                            bins_score[item][item2]+=1
                            bins_score_delta[item][item2]+=float(delta_coverage_vari)
                elif ave_cov > 15 and ave_cov <= 30:
                    if delta_coverage_vari <= 20: ### variable percentage
                        if item2 not in bins_score[item].keys():
                            bins_score[item][item2]=1
                            bins_score_delta[item][item2]=float(delta_coverage_vari) 
                        else:
                            bins_score[item][item2]+=1
                            bins_score_delta[item][item2]+=float(delta_coverage_vari)
                else:
                    if delta_coverage_vari <= 10: ### variable percentage
                        if item2 not in bins_score[item].keys():
                            bins_score[item][item2]=1
                            bins_score_delta[item][item2]=float(delta_coverage_vari) 
                        else:
                            bins_score[item][item2]+=1
                            bins_score_delta[item][item2]+=float(delta_coverage_vari)
            i+=1
        
        for item2 in bins_score[item].keys():
            if int(bins_score[item][item2]) < num:
                del bins_score[item][item2]
                del bins_score_delta[item][item2]

    binset_coverage_gc={}
    for item in binset_GC_ratio1.keys():
        for item2 in binset_GC_ratio2.keys():
            delta_GC_ratio=round(100*abs(float(binset_GC_ratio2[item2])-float(binset_GC_ratio1[item]))/float(binset_GC_ratio1[item]),2)
            # delta_GC_ratio=round(100*float(delta_GC)/float(binset_GC_ratio1[item]),1)
            if delta_GC_ratio <= 5: ## +- 3% total 6% var
                if item2 in bins_score[item].keys():
                    if str(item) not in binset_coverage_gc.keys():
                        binset_coverage_gc[str(item)]=str(item2)+':'+str(bins_score_delta[item][item2])+':'+str(delta_GC_ratio)
                    else:
                        binset_coverage_gc[str(item)]+='\t'+str(item2)+':'+str(bins_score_delta[item][item2])+':'+str(delta_GC_ratio)

                if item not in bin_gc.keys():
                    bin_gc[item]=str(item)+':'+str(binset_GC_ratio1[item])+'\t'+str(item2)+':'+str(binset_GC_ratio2[item2])+':'+str(delta_GC_ratio)
                else:
                    bin_gc[item]+='\t'+str(item2)+':'+str(binset_GC_ratio2[item2])+':'+str(delta_GC_ratio)

    bins_coverage_score={}
    f_bin_coverage=open('Bins_similar_coverage_'+str(iteration_num)+'.txt', 'w')
    f_bin_coverage.write('Target Bin'+'\t'+'Similar Bins'+'\t'+'Average Coverage Variation'+'\n')
    for item in bins_score.keys():
        for item2 in bins_score_delta[item].keys():
            average_coverage_var=round(float(bins_score_delta[item][item2]/num),2)
            f_bin_coverage.write(str(item)+'\t'+str(item2)+'\t'+str(average_coverage_var)+'\n')
            bins_coverage_score[str(item)+'\t'+str(item2)]=str(average_coverage_var)
    f_bin_coverage.close()

    f_gc=open('Bins_similar_GC_'+str(iteration_num)+'.txt', 'w')
    f_gc.write('Group'+'\t'+'Bin:GC%:deltaGC%'+'\n')
    n=0
    for item in bin_gc.keys():
        n+=1
        f_gc.write(str(n)+'\t'+str(bin_gc[item])+'\n')
    f_gc.close()

    f_co=open('Bins_similar_co_GC_coverage_'+str(iteration_num)+'.txt', 'w')
    f_co.write('Target Bin'+'\t'+'Similar:Coverage_var(%):GC_var(%)'+'\n')
    n=0
    for item in binset_coverage_gc.keys():
        n+=1
        f_co.write(str(item)+'\t'+str(binset_coverage_gc[item])+'\n')
    f_co.close()
    print '-----------------------'
    return bins_coverage_score, bin_gc, binset_coverage_gc

def seq_comparitor(blast_output, binset1, binset2, seq_len1, seq_len2, binset_record1, binset_record2, binset_genome_size1, binset_genome_size2, bins_coverage_score):
    print 'Comparing blast output of', binset1, 'and', binset2
    print '-----------------------'
    contig_contig_score={}
    contig_contig_precentage={}
    contig_bin_genome={}
    contig_contig_score_mock={}

    binset_record1_mock={}
    binset_record1_mock[str(binset1)]={}
    binset_record2_mock={}
    binset_record2_mock[str(binset2)]={}
    for item in binset_record1[str(binset1)].keys():
        binset_record1_mock[str(binset1)][item]=0
    for item in binset_record2[str(binset2)].keys():
        binset_record2_mock[str(binset2)][item]=0

    # print binset_record1
    print 'Parsing', blast_output
    for line in open(blast_output,'r'):
        simi=str(line).strip().split('\t')[2]
        if float(simi) >= 98.5:
            ID1_ORFs=str(line).strip().split('\t')[1] ### Be aware of the difference between ID1 and ID2. Seq1 has been set as database. 
            ID1=ID1_ORFs.split('_')[0]
            ID2_ORFs=str(line).strip().split('\t')[0]
            ID2=ID2_ORFs.split('_')[0]
            length=str(line).strip().split('\t')[3]
            Alignment_length=int(int(length)*float(simi)/100)
            ID1_aligned_percentage=round(float(int(length)*float(simi)/seq_len1[ID1_ORFs]),2)
            ID2_aligned_percentage=round(float(int(length)*float(simi)/seq_len2[ID2_ORFs]),2)
            num=len(contig_contig_score_mock)
            contig_contig_score_mock[str(ID1)+'\t'+str(ID2)]=1
            if len(contig_contig_score_mock) == num+1:
                contig_contig_score[str(ID1)+'\t'+str(ID2)]=Alignment_length
                contig_contig_precentage[str(ID1)+'\t'+str(ID2)]=str(ID1_aligned_percentage)+'\t'+str(ID2_aligned_percentage)
            else:
                contig_contig_score[str(ID1)+'\t'+str(ID2)]+=Alignment_length
                ID1_aligned_percentage_2=float(contig_contig_precentage[str(ID1)+'\t'+str(ID2)].split('\t')[0])+ID1_aligned_percentage
                ID2_aligned_percentage_2=float(contig_contig_precentage[str(ID1)+'\t'+str(ID2)].split('\t')[1])+ID2_aligned_percentage
                contig_contig_precentage[str(ID1)+'\t'+str(ID2)]=str(ID1_aligned_percentage_2)+'\t'+str(ID2_aligned_percentage_2)

            bin_name=[]
            m1=len(binset_record1_mock[str(binset1)])
            binset_record1_mock[str(binset1)][str(ID1)]=0
            if len(binset_record1_mock[str(binset1)]) == m1:
                m2=len(binset_record2_mock[str(binset2)])
                binset_record2_mock[str(binset2)][str(ID2)]=0
                if len(binset_record2_mock[str(binset2)]) == m2:
                    bin1=str(binset_record1[str(binset1)][str(ID1)])
                    bin2=str(binset_record2[str(binset2)][str(ID2)])
                    if ',' not in bin1 and ',' not in bin2:
                        bin_name.append(bin1+'\t'+bin2)
                    elif ',' in bin1 and ',' not in bin2:
                        lis=str(bin1).split(',')
                        xx=0
                        while xx < len(lis):
                            bin_name.append(str(lis[xx])+'\t'+bin2)
                            xx+=1
                    elif ',' not in bin1 and ',' in bin2:
                        lis=bin2.split(',')
                        xx=0
                        while xx < len(lis):
                            bin_name.append(bin1+'\t'+bin2.split(',')[xx])
                            xx+=1
                    else:
                        lis1=bin1.split(',')
                        lis2=bin2.split(',')
                        for item in lis1:
                            xx=0
                            while xx < len(lis2):
                                bin_name.append(str(item)+'\t'+bin2.split(',')[xx])
                                xx+=1
                    
                    for item in bin_name:
                        if str(item) not in contig_bin_genome.keys():
                            contig_bin_genome[str(item)]=Alignment_length
                        else:
                            contig_bin_genome[str(item)]+=Alignment_length
                else:
                    del binset_record2_mock[str(binset2)][str(ID2)]
            else:
                del binset_record1_mock[str(binset1)][str(ID1)]

    print 'Sorting contigs score'
    contig_contig_score=sorted(contig_contig_score.items(), key=lambda d:d[0])

    f=open('Contigs_scoring_'+blast_output, 'w')
    f.write('Contig1'+'\t'+'Contig2'+'\t'+'Alignment_length'+'\t'+'Contig1_aligned_percentage'+'\t'+'Contig2_aligned_percentage'+'\n')
    for item in contig_contig_score:
        f.write(str(item[0])+'\t'+str(item[1])+'\t'+str(contig_contig_precentage[item[0]])+'\n')
    f.close()

    for item in contig_bin_genome.keys():
        num=contig_bin_genome[item]
        A_percentage=round(100*float(num)/float(binset_genome_size1[str(binset1)][item.split('\t')[0]]),2)
        B_percentage=round(100*float(num)/float(binset_genome_size2[str(binset2)][item.split('\t')[1]]),2)
        contig_bin_genome[item]=str(num)+'\t'+str(A_percentage)+'\t'+str(B_percentage)

    contig_bin_genome=sorted(contig_bin_genome.items(), key=lambda d:d[0])

    bin_seq_coverage_filtrated={}
    f=open('Raw_Bins_Comparison_'+blast_output, 'w')
    # f2=open('Filtrated_Bins_Comparison_'+blast_output, 'w')
    f3=open('Filtrated_Bins_Seq_Similarity_Coverage_'+blast_output, 'w')
    title='GenomeA'+'\t'+'GnomesB'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'
    f.write(title+'\n')
    # f2.write(title+'\n')
    f3.write(title+'\t'+'Average_Coverage_Variation'+'\n')
    for item in contig_bin_genome:
        f.write(str(item[0])+'\t'+str(item[1])+'\n')
        # print str(item[1])
        alignment_length=int(str(item[1].split('\t')[0]))
        A_sim=item[1].split('\t')[1]
        B_sim=item[1].split('\t')[2]
        Sum_sim=float(A_sim)+float(B_sim)
        if 'noclass' not in item[0] and 'unbined' not in item[0] and alignment_length >= 500000 and Sum_sim >= 100: ### Similar sequence which is larger than 100 Kbp will be kept to consider
            if str(item[0]) in bins_coverage_score.keys():
                    # f2.write(str(item[0])+'\t'+str(item[1])+'\n')
                bin_seq_coverage_filtrated[str(item[0])]=str(item[1])+'\t'+str(bins_coverage_score[str(item[0])])
                f3.write(str(item[0])+'\t'+str(item[1])+'\t'+str(bins_coverage_score[str(item[0])])+'\n')
            else:
                continue
    f.close()
    # f2.close()
    f3.close()
    return contig_contig_score, contig_bin_genome, bin_seq_coverage_filtrated

def checkm_connections(binset):
    print 'Reading checkm output of', binset
    pwd=os.getcwd()
    binset_checkm_connection={}
    for root, dirs, files in os.walk(binset):
        os.chdir(pwd+'/'+binset)
        for file in files:
            if '_bin_stats_ext.tsv' in file:
                for line in open(file, 'r'):
                    bins=str(line).strip().split('\t')[0]
                    if '_maxbin2_genomes.' in bins or '_concoct_genomes.' in bins:
                        bin_id=str(bins)+'.fasta'
                    elif '_metabat_genomes.' in bins:
                        bin_id=str(bins)+'.fa'
                    else:
                        continue
                    binset_checkm_connection[str(bin_id)]={}
                    binset_checkm_connection[str(bin_id)]['Connections']=int(str(line).strip().split('Connections\': ')[1].split(',')[0])
                    binset_checkm_connection[str(bin_id)]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('\',')[0]
                    binset_checkm_connection[str(bin_id)]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split(',')[0])
                    binset_checkm_connection[str(bin_id)]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split(',')[0])
                    binset_checkm_connection[str(bin_id)]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0])
    os.chdir(pwd)
    print 'Done of reading checkm output of', binset
    print '-----------------------------------------'
    return binset_checkm_connection

def bin_comparitor(bin_seq_coverage_filtrated, binset_checkm_connection_1, binset_checkm_connection_2, num):
    print 'Comparing bins'
    selected_bins, eliminated_bins={}, {}
    marker_score={'root':0, 'k':1, 'p':1.5, 'c':2.3, 'o':3.4, 'f':5.1, 'g':7.6, 's':11.4}
    f=open('Selected_bins_'+str(num)+'.txt', 'w')
    f.write('Selected_bin'+'\t'+'Similar_bins'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'+'\t'+'Average_Coverage_Variation'+'\t'+'Bin1_factors'+'\t'+'Bin2_factors'+'\n')
    for item in bin_seq_coverage_filtrated.keys():
        ID1=str(item).split('\t')[0]
        ID2=str(item).split('\t')[1]
        bin_score_1=marker_score[binset_checkm_connection_1[ID1]['marker lineage'].split('__')[0]]
        bin_score_2=marker_score[binset_checkm_connection_2[ID2]['marker lineage'].split('__')[0]]
        CPN_CTN_1=binset_checkm_connection_1[ID1]['Completeness']-binset_checkm_connection_1[ID1]['Contamination']
        CPN_CTN_2=binset_checkm_connection_2[ID2]['Completeness']-binset_checkm_connection_2[ID2]['Contamination']
        if bin_score_1 == bin_score_2:    
            if CPN_CTN_1 == CPN_CTN_2:
                bin1=binset_checkm_connection_1[ID1]['Connections']/binset_checkm_connection_1[ID1]['Genome size']
                bin2=binset_checkm_connection_2[ID2]['Connections']/binset_checkm_connection_2[ID2]['Genome size']
                if bin1 <= bin2:
                # delta_connections_1_2=binset_checkm_connection_1[ID1]['Connections']-binset_checkm_connection_2[ID2]['Connections']
                # delta_genome_size_1_2=binset_checkm_connection_1[ID1]['Genome size']-binset_checkm_connection_2[ID2]['Genome size']
                # bin_delta=(delta_connections_1_2/delta_genome_size_1_2)-(binset_checkm_connection_2[ID2]['Connections']/binset_checkm_connection_2[ID2]['Genome size'])
                # if bin_delta >= 0:
                    selected_bins[ID1]=binset_checkm_connection_1[ID1]
                    eliminated_bins[ID2]=binset_checkm_connection_2[ID2]
                    f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
                else:
                    selected_bins[ID2]=binset_checkm_connection_2[ID2]
                    eliminated_bins[ID1]=binset_checkm_connection_1[ID1]
                    f.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
            elif CPN_CTN_1 > CPN_CTN_2:
                selected_bins[ID1]=binset_checkm_connection_1[ID1]
                eliminated_bins[ID2]=binset_checkm_connection_2[ID2]
                f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
            else:
                selected_bins[ID2]=binset_checkm_connection_2[ID2]
                eliminated_bins[ID1]=binset_checkm_connection_1[ID1]
                f.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
        else:
            if float(bin_score_1)*float(CPN_CTN_1) >= float(bin_score_2)*float(CPN_CTN_2):
                selected_bins[ID1]=binset_checkm_connection_1[ID1]
                eliminated_bins[ID2]=binset_checkm_connection_2[ID2]
                f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
            else:
                selected_bins[ID2]=binset_checkm_connection_2[ID2]   
                eliminated_bins[ID1]=binset_checkm_connection_1[ID1]
                f.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n') 
    f.close()

    print 'Comparison done!'
    print '----------------'
    return selected_bins, eliminated_bins

def new_selected_bins_generator(selected_bins, eliminated_bins, all_bins_1, all_bins_2, binset_checkm_connection_1, binset_checkm_connection_2, iteration_num, binset1, binset2, ORFs_assembly1, ORFs_assembly2, coverage_matrix1, coverage_matrix2):
    pwd=os.getcwd()
    extract_bins_1, extract_bins_2, coverage_maxtrix={}, {}, {}

    nx=0
    for line in open(coverage_matrix1, 'r'):
        nx+=1
        if nx == 1:
            coverage_title=str(line).strip()
        else:
            coverage_maxtrix[str(line).strip().split('\t')[0]]=str(line).strip()
    
    nx=0
    for line in open(coverage_matrix2, 'r'):
        nx+=1
        if nx > 1:
            coverage_maxtrix[str(line).strip().split('\t')[0]]=str(line).strip()

    for item in all_bins_1.keys():
        if item not in selected_bins.keys() and item not in eliminated_bins.keys():
            extract_bins_1[item]=binset_checkm_connection_1[item]
    
    for item in all_bins_2.keys():
        if item not in selected_bins.keys() and item not in eliminated_bins.keys():
            extract_bins_2[item]=binset_checkm_connection_2[item]

    f=open('Extract_bins_'+str(iteration_num)+'.txt','w')
    f.write('Extract bins from binset 1'+'\n')
    for item in extract_bins_1.keys():
        f.write(str(item)+'\t'+str(extract_bins_1[item])+'\n')

    f.write('Extract bins from binset 2'+'\n')
    for item in extract_bins_2.keys():
        f.write(str(item)+'\t'+str(extract_bins_2[item])+'\n')
    f.close()

    try:
        os.mkdir('Iteration_'+str(iteration_num)+'_genomes')
    except:
        os.system('rm -rf Iteration_'+str(iteration_num)+'_genomes')
        os.mkdir('Iteration_'+str(iteration_num)+'_genomes')
        print 'Iteration_'+str(iteration_num)+'_genomes exist'
        print 'Re-created folder of Iteration_'+str(iteration_num)+'_genomes'

    new_Contigs, new_Contig_mock={}, {}
    for root, dirs, files in os.walk(binset1):
        os.chdir(binset1)
        for file in files:
            if file in extract_bins_1.keys() or file in selected_bins.keys():
                os.system('cp '+file+' '+pwd+'/Iteration_'+str(iteration_num)+'_genomes')
                for record in SeqIO.parse(file, 'fasta'):
                    new_Contigs[record.id]=str(record.seq)
                    new_Contig_mock[record.id]=0

    for root, dirs, files in os.walk(pwd+'/'+binset2):
        os.chdir(pwd+'/'+binset2)
        for file in files:
            if file in extract_bins_2.keys() or file in selected_bins.keys():
                os.system('cp '+file+' '+pwd+'/Iteration_'+str(iteration_num)+'_genomes')
                for record in SeqIO.parse(file, 'fasta'):
                    new_Contigs[record.id]=str(record.seq)
                    new_Contig_mock[record.id]=0

    os.chdir(pwd)                    
    new_contigs_name='Contigs_iteration_'+str(iteration_num)+'.fa'
    new_contigs=open(new_contigs_name, 'w')
    new_coverage_name='Coverage_matrix_for_binning_iteration_'+str(iteration_num)+'.txt'
    new_coverage=open(new_coverage_name,'w')
    new_coverage.write(str(coverage_title)+'\n')
    
    for item in new_Contigs.keys():
        new_contigs.write('>'+str(item)+'\n'+str(new_Contigs[item])+'\n')
        new_coverage.write(str(coverage_maxtrix[str(item)])+'\n')
    new_coverage.close()
    new_contigs.close()

    # new_ORFs=ORFs_predictor(str(new_contigs_name))

    new_ORFs=open('ORFs_iteration_'+str(iteration_num)+'.fa','w')
    id_num=len(new_Contig_mock)
    for record in SeqIO.parse(ORFs_assembly1, 'fasta'):
        new_Contig_mock[str(record.id).split('_')[0]]=0
        if len(new_Contig_mock) == id_num:
        # if str(record.id).split('_')[0] in new_Contigs.keys():
           new_ORFs.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
        else:
            del new_Contig_mock[str(record.id).split('_')[0]]
    
    for record in SeqIO.parse(ORFs_assembly2, 'fasta'):
        new_Contig_mock[str(record.id).split('_')[0]]=0
        if len(new_Contig_mock) == id_num:
        # if str(record.id).split('_')[0] in new_Contigs.keys():
           new_ORFs.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
        else:
            del new_Contig_mock[str(record.id).split('_')[0]]
    new_ORFs.close()
    
    
    high_quality_bins,all_remain_bins={},{}
    os.chdir(pwd+'/Iteration_'+str(iteration_num)+'_genomes')
    f=open('Iteration_'+str(iteration_num)+'_bin_stats_ext.tsv','w')
    for item in selected_bins.keys():
        all_remain_bins[str(item)]=str(selected_bins[item])
        #if '_metabat_genomes.' in str(item):
        #    bin_name=str(item).replace('.fa', '')
        #elif '_maxbin2_genomes.' in str(item):
        #    bin_name=str(item).replace('.fasta', '')  
        #else:
        #    continue
        item_list=item.split('.')
        item_list.remove(item_list[-1])
        bin_name='.'.join(item_list)
        f.write(str(bin_name)+'\t'+str(selected_bins[item])+'\n')

    for item in extract_bins_1.keys():
        all_remain_bins[str(item)]=str(extract_bins_1[item])
        item_list=item.split('.')
        item_list.remove(item_list[-1])
        bin_name='.'.join(item_list)
        f.write(str(bin_name)+'\t'+str(extract_bins_1[item])+'\n')

    for item in extract_bins_2.keys():
        all_remain_bins[str(item)]=str(extract_bins_2[item])
        item_list=item.split('.')
        item_list.remove(item_list[-1])
        bin_name='.'.join(item_list)
        f.write(str(bin_name)+'\t'+str(extract_bins_2[item])+'\n')
    f.close()

    f=open('Medium_quality_bins_CPN-5CTN_50_iteration_'+str(iteration_num)+'_.txt','w')
    f2=open('High_quality_bins_CPN-5CTN_70_iteration_'+str(iteration_num)+'_.txt','w')
    for item in all_remain_bins.keys():
        cpn=str(all_remain_bins[item]).split('Completeness\': ')[1].split(',')[0]
        # if str(all_remain_bins[item]).count(',') == 4:
        ctn=str(all_remain_bins[item]).split('Contamination\': ')[1].split('}')[0]
        quality=float(cpn)-5*float(ctn)
        if quality >= 50:
            f.write(str(item)+'\t'+str(all_remain_bins[str(item)])+'\n')
        else:
            continue
        
        if quality >= 70:
            f2.write(str(item)+'\t'+str(all_remain_bins[str(item)])+'\n')
    f.close()
    f2.close()
    os.chdir(pwd)
    return extract_bins_1, extract_bins_2, 'ORFs_iteration_'+str(iteration_num)+'.fa', 'Iteration_'+str(iteration_num)+'_genomes', new_coverage_name

def multiple_assembly_comparitor(ORFs_assembly1, ORFs_assembly2, assembly1_binset, assembly2_binset, coverage_matrix1, coverage_matrix2, alignment_output, num): 
    binset_dict={}
    seq_len1=ORFs_sequence_length_recorder(ORFs_assembly1)
    seq_len2=ORFs_sequence_length_recorder(ORFs_assembly2)

    binset_record, binset_genome_size={}, {}
    A_group=genome_contigs_recorder(assembly1_binset, binset_record, binset_genome_size, coverage_matrix1)
    binset_record1=A_group[0]
    binset_genome_size1=A_group[1]
    binset_coverage_avg1=A_group[2]
    binset_GC_ratio1=A_group[3]
    all_bins_1=A_group[4]

    B_group=genome_contigs_recorder(assembly2_binset, binset_record, binset_genome_size, coverage_matrix2)
    binset_record2=B_group[0]
    binset_genome_size2=B_group[1]
    binset_coverage_avg2=B_group[2]
    binset_GC_ratio2=B_group[3]
    all_bins_2=B_group[4]

    C=coverage_GC_comparitor(binset_coverage_avg1, binset_coverage_avg2, binset_GC_ratio1, binset_GC_ratio2, num)
    bins_coverage_score=C[0]
    
    D=seq_comparitor(alignment_output, assembly1_binset, assembly2_binset, seq_len1, seq_len2, binset_record1, binset_record2, binset_genome_size1, binset_genome_size2, bins_coverage_score)
    bin_seq_coverage_filtrated=D[2]

    binset_checkm_connection_1=checkm_connections(assembly1_binset)
    binset_checkm_connection_2=checkm_connections(assembly2_binset)

    E=bin_comparitor(bin_seq_coverage_filtrated, binset_checkm_connection_1, binset_checkm_connection_2, num)
    selected_bins=E[0]
    eliminated_bins=E[1]

    F=new_selected_bins_generator(selected_bins, eliminated_bins, all_bins_1, all_bins_2, binset_checkm_connection_1, binset_checkm_connection_2, num, assembly1_binset, assembly2_binset, ORFs_assembly1, ORFs_assembly2, coverage_matrix1, coverage_matrix2)    
    new_ORFs=F[2]
    new_bin_folder=F[3]
    new_coverage=F[4]
    return new_ORFs, new_bin_folder, new_coverage

def final_binset_comparitor(ORFs, final_iteration_folder):
    pwd=os.getcwd()
    nr_ORFs, nr_ORFs_mock, ORFs_bin, bins={},{}, {}, {}
    os.chdir(pwd+'/'+final_iteration_folder)
    print 'Parsing bins'
    for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
        for file in files:
            if '.fasta' in file or '.fa' in file:
                bins[file]={}
                bins[file]['totallen']=0
                bins[file]['seq']={}
                for record in SeqIO.parse(file, 'fasta'):
                    ORFs_bin[record.id]=file
                    bins[file]['seq'][record.id]=len(record.seq)
                    bins[file]['totallen']+=len(record.seq)

    os.chdir(pwd)
    print 'Finding redanbunce ORFs'
    for record in SeqIO.parse(ORFs, 'fasta'):
        num=len(nr_ORFs_mock)
        nr_ORFs_mock[record.seq]=1
        if len(nr_ORFs_mock) > num:
            if '*' in str(record.seq):
                nr_ORFs[str(record.seq).replace('*','').strip()]=[ORFs_bin[str(record.id).split('_')[0]]]
            else:
                nr_ORFs[str(record.seq)]=[ORFs_bin[str(record.id).split('_')[0]]]
        else:
            if '*' in str(record.seq):
                nr_ORFs[str(record.seq).replace('*','').strip()].append(ORFs_bin[str(record.id).split('_')[0]])
            else:
                nr_ORFs[record.seq].append(ORFs_bin[str(record.id).split('_')[0]])

    f_re=open('Contigs_sharing_bins.txt', 'w')
    f_re.write('Sequence'+'\t'+'Bins'+'\n')
    for item in nr_ORFs.keys():
        if len(nr_ORFs[item]) == 1:
            del nr_ORFs[item]
        else:
            f_re.write(str(item)+'\t'+str(nr_ORFs[item])+'\n')
            # continue
    f_re.close()

    bin_similar={}
    print 'Finding bins sharing contigs'
    for item in nr_ORFs.keys():
        num=len(nr_ORFs[item])
        nr_ORFs[item].sort()
        for i in range(0, num-1):
            for k in range(i+1, num):
                item2=str(nr_ORFs[item][i])+'\t'+str(nr_ORFs[item][k])
                if item2 not in bin_similar.keys():
                    bin_similar[item2]=len(item)
                else:
                    bin_similar[item2]+=len(item)

    bin_checkm, bin_checkm_o={}, {}
    os.chdir(pwd+'/'+final_iteration_folder)
    for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
        for file in files:
            if '_bin_stats_ext.tsv' in file:
                for line in open(file, 'r'):
                    bin_id=str(line).strip().split('\t')[0]
                    if '_maxbin2_genomes.' in bin_id or '_concoct_genomes.' in bin_id:
                        bin_id_f=bin_id+'.fasta'
                    elif '_metabat_genomes.' in bin_id:
                        bin_id_f=bin_id+'.fa'
                    else:
                        continue
                    bin_checkm_o[bin_id]=str(line).strip()
                    bin_checkm[bin_id_f]={}
                    bin_checkm[bin_id_f]['Connections']=int(str(line).strip().split('Connections\': ')[1].split(',')[0])
                    bin_checkm[bin_id_f]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('\',')[0]
                    bin_checkm[bin_id_f]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split(',')[0])
                    bin_checkm[bin_id_f]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split(',')[0])
                    bin_checkm[bin_id_f]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0])

    os.chdir(pwd)
    try:
        os.system('mkdir BestBinset')
    except:
        os.system('rm -rf BestBinset')
        os.system('mkdir BestBinset')
        print 'BestBinset exist'
        print 'Re-created folder of BestBinset'
    
    os.chdir('BestBinset')

    filtrated_bin={}
    f=open('Sim_bin_in_final_iteration.txt','w')
    f1=open('Filtrated_high_sim_bin_in_final_iteration.txt','w')
    f.write('Bin1'+'\t'+'Bin2'+'\t'+'Aligned(bp)'+'\t'+'Percetage_bin1(%)'+'\t'+'Percetage_bin2(%)'+'\n')
    f1.write('Bin1'+'\t'+'Bin2'+'\t'+'Aligned(bp)'+'\t'+'Percetage_bin1(%)'+'\t'+'Percetage_bin2(%)'+'\n')
    for item in bin_similar.keys():
        bin1=item.split('\t')[0]
        bin2=item.split('\t')[1]
        pb1=round(100*int(bin_similar[item])/int(bins[bin1]['totallen']), 2) # be aware of different between aa and nt
        pb2=round(100*int(bin_similar[item])/int(bins[bin2]['totallen']), 2)
        f.write(item+'\t'+str(bin_similar[item])+'\t'+str(pb1)+'\t'+str(pb2)+'\n')
        sum_sim=pb1+pb2
        if sum_sim >= 140:
            f1.write(item+'\t'+str(bin_similar[item])+'\t'+str(pb1)+'\t'+str(pb2)+'\n')
            filtrated_bin[item]=bin_similar[item]
        elif pb1 >= 80 :
            f1.write(item+'\t'+str(bin_similar[item])+'\t'+str(pb1)+'\t'+str(pb2)+'\n')
            filtrated_bin[item]=bin_similar[item]
        elif pb2 >= 80:
            f1.write(item+'\t'+str(bin_similar[item])+'\t'+str(pb1)+'\t'+str(pb2)+'\n')
            filtrated_bin[item]=bin_similar[item]
    f.close()
    f1.close()
    
    selected_bins, eliminated_bins={}, {}
    marker_score={'root':0, 'k':1, 'p':2, 'c':3, 'o':4, 'f':5, 'g':6, 's':7}
    f=open('Selected_bins_best_binset.txt', 'w')
    f.write('Selected_bin'+'\t'+'Similar_bins'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'+'\t'+'Average_Coverage_Variation'+'\t'+'Bin1_factors'+'\t'+'Bin2_factors'+'\n')
    for item in filtrated_bin.keys():
        ID1=str(item).split('\t')[0]
        ID2=str(item).split('\t')[1]
        bin_score_1=marker_score[bin_checkm[ID1]['marker lineage'].split('__')[0]]
        bin_score_2=marker_score[bin_checkm[ID2]['marker lineage'].split('__')[0]]
        if bin_score_1 == bin_score_2:
            CPN_CTN_1=bin_checkm[ID1]['Completeness']-bin_checkm[ID1]['Contamination']
            CPN_CTN_2=bin_checkm[ID2]['Completeness']-bin_checkm[ID2]['Contamination']
            if CPN_CTN_1 == CPN_CTN_2:
                bin1=bin_checkm[ID1]['Connections']/bin_checkm[ID1]['Genome size']
                bin2=bin_checkm[ID2]['Connections']/bin_checkm[ID2]['Genome size']
                if bin1 <= bin2:
                # delta_connections_1_2=binset_checkm_connection_1[ID1]['Connections']-binset_checkm_connection_2[ID2]['Connections']
                # delta_genome_size_1_2=binset_checkm_connection_1[ID1]['Genome size']-binset_checkm_connection_2[ID2]['Genome size']
                # bin_delta=(delta_connections_1_2/delta_genome_size_1_2)-(binset_checkm_connection_2[ID2]['Connections']/binset_checkm_connection_2[ID2]['Genome size'])
                # if bin_delta >= 0:
                    selected_bins[ID1]=bin_checkm[ID1]
                    eliminated_bins[ID2]=bin_checkm[ID2]
                    f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                else:
                    selected_bins[ID2]=bin_checkm[ID2]
                    eliminated_bins[ID1]=bin_checkm[ID1]
                    f.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
            elif CPN_CTN_1 > CPN_CTN_2:
                selected_bins[ID1]=bin_checkm[ID1]
                eliminated_bins[ID2]=bin_checkm[ID2]
                f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
            else:
                selected_bins[ID2]=bin_checkm[ID2]
                eliminated_bins[ID1]=bin_checkm[ID1]
                f.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
        elif bin_score_1 > bin_score_2:
            selected_bins[ID1]=bin_checkm[ID1]
            eliminated_bins[ID2]=bin_checkm[ID2]
            f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
        else:
            selected_bins[ID2]=bin_checkm[ID2]   
            eliminated_bins[ID1]=bin_checkm[ID1]
            f.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n') 
    f.close()

    # for item in eliminated_bins.keys():
    #    if item in selected_bins.keys():
    #        del selected_bins[item]

    f=open('Best_binset_bin_stats_ext.tsv', 'w')
    f2=open('Possible_same_bin.txt','w')
    f3=open('Highly_possible_same_bin.txt','w')
    for item in bin_checkm_o.keys():
        if item not in eliminated_bins.keys():
            f.write(str(bin_checkm_o[item])+'\n')
            Completeness=str(bin_checkm_o[item]).split('Completeness\': ')[1].split(',')[0]
            Contamination=str(bin_checkm_o[item]).split('Contamination\': ')[1].split('}')[0]
            delta_value=float(Completeness)-5*float(Contamination)
            if delta_value >= 70:
                f3.write(str(bin_checkm_o[item])+'\n')
            elif delta_value >= 50:
                f2.write(str(bin_checkm_o[item])+'\n')
            else:
                continue
    f.close()
    f2.close()
    f3.close()

    f=open('Genome_group_all_list_BestBinsSet.txt','w')
    # f2=open('prebinned_genomes_output_for_dataframe_BestBinsSet.txt','w')
    # f3=open('BestBinsSet.depth.txt','w')

    print 'Parsing files in Best Binset'
    os.chdir(pwd+'/'+final_iteration_folder)
    for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('_genomes.')[1]
                if '.fasta' in hz or '.fa' in hz:
                    if file not in eliminated_bins.keys():
                        os.system('cp '+file+' '+pwd+'/BestBinset')
    
    # contig_bin={}
    os.chdir(pwd+'/BestBinset')
    # for root, dirs, files in os.walk(pwd+'/BestBinset'):
    #     for file in files:
    #         if '_genomes.' in file:
    #             hz=str(file).split('_genomes.')[1]
    #             if '.fasta' in hz or '.fa' in hz:
    #                 for record in SeqIO.parse(file, 'fasta'):
    #                     contig_bin[record.id]=file

    # title1, title2, title3=[], [], []
    # os.chdir(pwd+'/'+final_iteration_folder)
    # for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
    bin_folders, title={}, []
    for root, dirs, files in os.walk(pwd+'/BestBinset'):
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('_genomes.')[1]
                if '.fasta' in hz or '.fa' in hz:
                    bin_source_folder=str(file).split('_genomes.')[0]
                    bin_folders[bin_source_folder]=1
            # if '.depth.txt' in file:
            #     n=0
            #     for line in open(file, 'r'):
            #         n+=1
            #         if n == 1 and len(title1) == 0:
            #             title1.append(str(line))
            #             f3.write(str(line))
            #         else:
            #             if str(line).strip().split('\t')[0] in contig_bin.keys():
            #                 f3.write(str(line))

    for item in bin_folders.keys():
        os.chdir(pwd+'/'+item+'_genomes')
        for root, dirs, files in os.walk(pwd+'/'+item+'_genomes'):
            for file in files:
                if 'Genome_group_all_list_' in file:
                    n=0
                    for line in open(file, 'r'): 
                        n+=1
                        if n == 1 and len(title) == 0:
                            f.write(str(line))
                            title.append(str(line))
                        else:
                            if str(line).strip().split('\t')[0] not in eliminated_bins.keys():
                                f.write(str(line))
    f.close()

    os.chdir(pwd)
    os.system('mkdir BestBinset_comparison_files')
    os.system('rm temp.orfs.* Contigs_iteration_*')
    os.system('mv Contigs_sharing_bins.txt Contigs_iteration_* Coverage_matrix_for_binning_iteration_* Selected_bins_* Extract_bins_* Bins_similar_* Contigs_scoring_* Raw_bins_* Filtrated_Bins_* *_vs_* *_BestBinSet_gc.txt *_BestBinSet_bins_coverage.txt BestBinset_comparison_files')
            # if 'prebinned_genomes_output_for_dataframe_' in file:
            #     n=0
            #     for line in open(file, 'r'):
            #         n+=1
            #         if n == 1 and len(title2) == 0:
            #             title2.append(str(line))
            #             f2.write(str(line))
            #         else:
            #             if str(line).strip().split('\t')[1] in bin_selected.keys():
            #                 f2.write(str(line))
    
    # f2.close()
    # f3.close()          
    print 'Done!'

def multiple_assembly_comparitor_main(Contig_list, BestBinSet_list, Coverage_list):
    print 'Processing', Contig_list
    ORF_list=[]
    contig_num=len(Contig_list)
    for item in Contig_list:
       ORF_list.append(ORFs_predictor(item))
    #ORF_list=['1_9groups_assembly.fa.orfs.fna', '2_AN1-32_assembly.fa.orfs.fna', '3_AN1-52_assembly.fa.orfs.fna', '4_AN1-67D_assembly.fa.orfs.fna', '5_AN1-118_assembly.fa.orfs.fna', '6_AN1-146_assembly.fa.orfs.fna', '7_AN1-188_assembly.fa.orfs.fna', '8_AN1-167_assembly.fa.orfs.fna', '9_AN2-66_assembly.fa.orfs.fna', '10_AN2-109_assembly.fa.orfs.fna', '11_AN2-284-pacbio-contigs.fasta.orfs.fna']

    print 'Processing', BestBinSet_list

    print 'Processing', Coverage_list

    alignment_list=[]
    num=0
    for i in range(1, contig_num):
        alignment_output=ORFs_aligner(ORF_list[0], ORF_list[-1])
        alignment_list.append(alignment_output)
        new=multiple_assembly_comparitor(ORF_list[0], ORF_list[-1], BestBinSet_list[0], BestBinSet_list[-1], Coverage_list[0], Coverage_list[-1], alignment_output, i)
        ORF_list.remove(ORF_list[0])
        ORF_list.remove(ORF_list[-1])
        ORF_list.append(new[0])
        print 'Adding', str(new[0]), 'to ORFs list'
        BestBinSet_list.remove(BestBinSet_list[0])
        BestBinSet_list.remove(BestBinSet_list[-1])
        BestBinSet_list.append(new[1])
        print 'Adding', str(new[1]), 'to BestBinSet list'
        Coverage_list.remove(Coverage_list[0])
        Coverage_list.remove(Coverage_list[-1])
        Coverage_list.append(new[2])
        print 'Adding', str(new[2]), 'to coverage list'
    
    final_binset_comparitor(new[0], new[1])

if __name__ == '__main__': 
    Contig_list=['1_ZH_FAST_final.contigs.fa','2_ZH_MO_FAST_cat_final.contigs.fa', '3_ZH_MO_final.contigs.fa'] ### Beaware the ID of assembly should be different
    BestBinSet_list=['1_ZH_FAST_final.contigs.fa_BestBinsSet','2_ZH_MO_FAST_cat_final.contigs.fa_BestBinsSet', '3_ZH_MO_final.contigs.fa_BestBinsSet']
    Coverage_list=['Coverage_matrix_for_binning_1_ZH_FAST_final.contigs.fa.txt', 'Coverage_matrix_for_binning_2_ZH_MO_FAST_cat_final.contigs.fa.txt', 'Coverage_matrix_for_binning_3_ZH_MO_final.contigs.fa.txt']
    
    multiple_assembly_comparitor_main(Contig_list, BestBinSet_list, Coverage_list)
