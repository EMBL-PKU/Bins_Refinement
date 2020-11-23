#!/usr/bin/python
# -*- coding: UTF-8 -*-
#coding=utf-8

import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
from Bio import SeqIO
import os
from sklearn.decomposition import PCA
import numpy as np
# from Outlier_remover import *

def record_bin_coverage(binset):
    pwd=os.getcwd()
    bin_contigs, bin_contigs_mock, total_bin_contigs, assembly_list, m={}, {}, {}, {}, 0
    print 'Parsing bins'
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        os.chdir(pwd+'/'+str(binset))
        for file in files:
            if '_genomes.' in file:
                hz=file.split('_genomes.')[1]
                qz=file.split('_genomes.')[0]
                assembly_name_list=qz.split('_')
                assembly_name_list.remove(assembly_name_list[-1])
                assembly_name_list.remove(assembly_name_list[-1])
                assembly_name='_'.join(assembly_name_list)
                assembly_list[assembly_name]=1
                if '.fasta' in hz or '.fa' in hz:
                    m+=1
                    bin_contigs[file]={}
                    bin_contigs_mock[file]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_contigs[file][str(record.id)]=str(record.seq)
                        bin_contigs_mock[file][str(record.id)]=0
    print 'Parsed', m, 'bins'
    
    os.chdir(pwd)
    for assemblies in assembly_list.keys():
        for record in SeqIO.parse(assemblies, 'fasta'):
            total_bin_contigs[str(record.id)]=str(record.seq)

    print 'Recording the coverage of contigs from bins'
    contig_cov={}
    for assembly in assembly_list.keys():
        n=0
        for line in open('Coverage_matrix_for_binning_'+str(assembly)+'.txt','r'):
            n+=1
            if n >= 2:
                ls=str(line).strip().split('\t')
                num=(len(ls)-4)/3
                ids=str(line).strip().split('\t')[0]
                contig_cov[ids]={}
                for i in range(1,num+1):
                    contig_cov[ids][i]=float(str(line).strip().split('\t')[3*i+1])
        
    bin_contig_cov={}    
    for contig_id in contig_cov.keys():
        for bin_name in bin_contigs.keys():
            if bin_name not in bin_contig_cov.keys():
                bin_contig_cov[bin_name]={} 
            bin_contigs_mock[bin_name][contig_id]=1
            # if contig_id in bin_contigs[bin_name].keys():
            if len(bin_contigs_mock[bin_name]) == len(bin_contigs[bin_name]):
                bin_contig_cov[bin_name][contig_id]={} 
                for i in range(1,num+1):
                    bin_contig_cov[bin_name][contig_id][i]=contig_cov[contig_id][i]
            else:
                del bin_contigs_mock[bin_name][contig_id]
        
    print 'Writing bin coverage matrix file'
    os.system('mkdir Bin_coverage_after_contamination_removal')
    os.chdir('Bin_coverage_after_contamination_removal')
    for bins in bin_contig_cov.keys():
        f=open(bins+'_coverage_matrix.txt', 'w')
        f.write('Contig'+'\t'+'Coverage'+'\n')
        for contigs in bin_contig_cov[bins].keys():
            f.write(str(contigs)+'\t'+str(bin_contig_cov[bins][contigs])+'\n')
        f.close()
    os.chdir(pwd)
    return bin_contig_cov, bin_contigs, contig_cov, total_bin_contigs

def PE_connecting_contigs(assembly, PE_connections_file, binset, kmer_file):
    pwd=os.getcwd()
    print 'Finding common contigs of---', assembly, '---from---', binset, 'using ---', PE_connections_file
    bin_contigs, bin_contigs_mock, bin_seq, bin_select_contigs={}, {}, {}, {}
    os.chdir(pwd+'/'+binset)
    for root, dirs, files in os.walk(pwd+'/'+binset):
        # os.chdir(pwd+'/'+str(binset_1_assembly))
        for file in files:
            if '_genomes.' in file and str(assembly) in file:
                hz=file.split('_genomes.')[1]
                qz=file.split('_genomes.')[0]
                if '.fa' in hz or '.fasta' in hz:
                    # print 'Parsing', file
                    bin_contigs[file]={}
                    bin_select_contigs[file]={}
                    bin_contigs_mock[file]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_contigs[file][str(record.id)]=0
                        bin_select_contigs[file][str(record.id)]=0
                        bin_contigs_mock[file][str(record.id)]=0
                        bin_seq[str(record.id)]=str(record.seq)

    os.chdir(pwd)
    print 'Parsing', assembly, 'PE connections file'
    n, connections=0, {}
    for line in open(PE_connections_file,'r'):
        n+=1
        if n >= 2:
            id1=str(line).strip().split('\t')[0]
            id2=str(line).strip().split('\t')[2]
            connection_num=str(line).strip().split('\t')[3]
            if id1 not in connections.keys():
                connections[id1]={}
                connections[id1][str(id2)]=str(connection_num) ### Number suggests the connecting level
            else:
                connections[id1][str(id2)]=str(connection_num)

            if id2 not in connections.keys():
                connections[id2]={}
                connections[id2][str(id1)]=str(connection_num)
            else:
                connections[id2][str(id1)]=str(connection_num)

    print 'Parsing', assembly, 'kmer file'
    n, contig_bin_kmer, contig_kmer=0, {}, {}
    for line in open(kmer_file,'r'):
        n+=1
        if n >= 2:
            ids=str(line).strip().split('\t')[0]
            kmer_list=str(line).strip().split('\t')
            kmer_list.remove(kmer_list[0])
            kmer='\t'.join(kmer_list)
            contig_kmer[ids]=kmer
            for bin in bin_contigs.keys():
                if ids in bin_contigs[bin].keys():
                    if bin not in contig_bin_kmer.keys():
                        contig_bin_kmer[bin]={}
                        contig_bin_kmer[bin][ids]=kmer
                    else:
                        contig_bin_kmer[bin][ids]=kmer
    
    f=open('Connections_'+str(assembly)+'.txt', 'w')
    for item in connections.keys():
        f.write(str(item)+'\t'+str(connections[item])+'\n')
    f.close()

    bin_connecting_contigs, bin_connecting_contigs2={}, {}
    for bins in bin_contigs.keys():
        bin_connecting_contigs[bins]={}
        bin_connecting_contigs2[bins]={}
        bin_connecting_contigs[bins]=bin_contigs[bins]
        m, m_before, m_after=1, 0, 1
        while m <= 3: ### Consider connection level less than 3
            print bins, 'cycle', m
            if m_before != m_after:
                m_before = len(bin_connecting_contigs[bins])
                for contigs in connections.keys():
                    for item in connections[contigs].keys():
                        if contigs in bin_connecting_contigs[bins].keys() and item not in bin_connecting_contigs[bins].keys():
                            bin_connecting_contigs[bins][item]=m
                            bin_connecting_contigs2[bins][item]=m
                m_after = len(bin_connecting_contigs[bins])
                m+=1
            else:
                m=6

    for item in bin_connecting_contigs2.keys():
        if len(bin_connecting_contigs2[item]) == 0:
            del bin_connecting_contigs2[item]

    f=open('Bin_connecting_contigs_'+str(assembly)+'.txt', 'w')
    f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
    for item in bin_connecting_contigs2.keys():
        f.write(str(item)+'\t'+str(bin_connecting_contigs2[item])+'\n')
    f.close()

    try:
        os.mkdir('Bin_kmer')
    except:
        print 'Bin_kmer exists'

    os.chdir('Bin_kmer')
    for item in contig_bin_kmer.keys():
        f=open('Kmer_'+item, 'w')
        for contigs in contig_bin_kmer[item].keys():
            f.write(str(contigs)+'\t'+str(contig_bin_kmer[item][contigs])+'\n')
        f.close()
    os.chdir(pwd)

    return bin_connecting_contigs2, contig_bin_kmer, contig_kmer, connections

def test_outlier(connecting_contig, item_data, test_index):
    # print 'Judging', connecting_contig
    four = pd.Series(item_data).describe()
    # print(four)
    # print('Q1= {0}, Q2= {1}, Q3={2}'.format(four['25%'],four['50%'],four['75%']))
    Q1 = four['25%']
    Q3 = four['75%']
    IQR = Q3 - Q1
    upper1 = Q3 + 1.5 * IQR
    lower1 = Q1 - 1.5 * IQR
    stat=0
    if item_data[test_index] > float(upper1) or item_data[0] < float(lower1):
        stat=0
    else:
        stat=1
    return stat

def PCA_slector(data_array, num_contig):
    pca = PCA(n_components=1)
    pca.fit(data_array)
    explained_variance_ratio=pca.explained_variance_ratio_
    # print(explained_variance_ratio)
    # num=len(explained_variance_ratio)
    newData=pca.fit_transform(data_array)
    newData2=newData.reshape((1,num_contig))
    # print 'Shape', num, num_contig
    newData_list=newData2.tolist()
    n=0
    for item in newData_list:
        n+=1
        if n == 1:
            newData_list_item=item
    return newData_list_item, explained_variance_ratio

def checkm(bin_folder):
    pwd=os.getcwd()
    os.system('checkm lineage_wf -t 42 -x fa '+str(bin_folder)+' '+str(bin_folder)+'_checkm')
    os.chdir(str(bin_folder)+'_checkm/storage/')
    print 'Parsing '+bin_folder+' checkm output'
    refined_checkm={}
    for line in open('bin_stats_ext.tsv','r'):
        binID=str(line).strip().split('{\'')[0].strip()
        genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
        taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
        completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
        contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].strip()
        # GC=round(float(str(line).strip().split('\'GC\':')[1].split(', \'GCN4\'')[0].strip())*100, 1)
        refined_checkm[str(binID)]={}
        refined_checkm[str(binID)]['marker lineage']=str(taxon)
        refined_checkm[str(binID)]['Completeness']=float(completeness)
        refined_checkm[str(binID)]['Genome size']=float(genome_size)
        refined_checkm[str(binID)]['Contamination']=float(contamination)
    os.chdir(pwd)
    return refined_checkm

def bin_comparison(original_bin_folder, new_bins_checkm, new_bin_folder):
    pwd=os.getcwd()
    print 'Comparing bins before and after refining process'
    os.chdir(pwd+'/'+str(original_bin_folder))
    bin_checkm={}
    for root, dirs, files in os.walk(pwd+'/'+str(original_bin_folder)):
        for file in files:
            if '_bin_stats_ext.tsv' in file:
                for line in open(file, 'r'):
                    binID=str(line).strip().split('{\'')[0].strip()
                    # if binID+'.recruited' in refined_checkm.keys():
                    connections=int(str(line).strip().split(':')[1].split(',')[0].strip())
                    taxon=str(line).strip().split(': \'')[1].split('\'')[0].strip()
                    completeness=float(str(line).strip().split('Completeness\':')[1].split(',')[0].strip())
                    contamination=float(str(line).strip().split('Contamination\':')[1].split('}')[0].strip())
                    genome_size=float(str(line).strip().split('Genome size\':')[1].split(',')[0].strip())
                    bin_checkm[str(binID)]={}
                    bin_checkm[str(binID)]['Connections']=int(connections)
                    bin_checkm[str(binID)]['marker lineage']=str(taxon)
                    bin_checkm[str(binID)]['Completeness']=float(completeness)
                    bin_checkm[str(binID)]['Genome size']=float(genome_size)
                    bin_checkm[str(binID)]['Contamination']=float(contamination)

    os.chdir(pwd+'/'+str(new_bin_folder))
    bin_comparison, bin_comparison_list={}, {}
    for refined_bin in new_bins_checkm.keys():
        bin_id=refined_bin.split('.retrieved')[0]
        if bin_id in bin_checkm.keys():
            if bin_id not in bin_comparison.keys():
                bin_comparison[bin_id]=bin_id+'---'+refined_bin
                bin_comparison_list[bin_id]=[]
                bin_comparison_list[bin_id].append(bin_id)
                bin_comparison_list[bin_id].append(refined_bin)
            else:
                bin_comparison[bin_id]+='---'+refined_bin
                bin_comparison_list[bin_id].append(refined_bin)

    print bin_comparison_list
    f=open('Refined_bin_comparison.txt','w')
    f.write('Bin'+'\t'+'Related_bin'+'\t'+'checkm'+'\n')
    bestbin={}
    for item in bin_comparison_list.keys():
        bestbin[item]={}
        write_out=str(item)+'\t'+(bin_comparison[item])
        for i in range(0, len(bin_comparison_list[item])):
            if i == 0:
                bin_id=bin_comparison_list[item][i]
                write_out+='\t'+str(bin_checkm[bin_id])
                bestbin[item][bin_id]={}
                bestbin[item][bin_id]=bin_checkm[bin_id]
            else:
                refined_bin_id=bin_comparison_list[item][i]
                write_out+='---'+str(new_bins_checkm[refined_bin_id])
                re_connections=int(new_bins_checkm[refined_bin_id]['Connections'])
                re_cpn=float(new_bins_checkm[refined_bin_id]['Completeness'])
                re_ctn=float(new_bins_checkm[refined_bin_id]['Contamination'])
                re_taxon=str(new_bins_checkm[refined_bin_id]['marker lineage'])
                re_genome_size=float(new_bins_checkm[refined_bin_id]['Genome size'])
                re_delta=re_cpn-re_ctn
                for bin_id2 in bestbin[item].keys():
                    ori_bin=bin_id2 
                    ori_connections=int(bestbin[item][bin_id2]['Connections'])
                    ori_cpn=float(bestbin[item][bin_id2]['Completeness'])
                    ori_ctn=float(bestbin[item][bin_id2]['Contamination'])
                    ori_taxon=str(bestbin[item][bin_id2]['marker lineage'])
                    ori_genome_size=float(bestbin[item][bin_id2]['Genome size'])
                    ori_delta=ori_cpn-ori_ctn
                
                if re_delta > ori_delta:
                    del bestbin[item][ori_bin]
                    bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]
                elif re_delta == ori_delta:
                    if re_connections < ori_connections:
                        del bestbin[item][ori_bin]
                        bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]
                else:
                    continue
        f.write(str(write_out)+'\n')
    f.close()

    f=open('Selected_bin.txt', 'w')
    f.write('Original bin'+'\t'+'Selected bin'+'\t'+'Factors'+'\n')
    selected_bin={}
    for item in bestbin.keys():
        for bins in bestbin[item].keys():
            f.write(str(item)+'\t'+str(bins)+'\t'+str(bestbin[item][bins])+'\n')
            selected_bin[bins]=1
    f.close()

    n=0
    for root, dirs, files in os.walk(pwd+'/'+str(new_bin_folder)):
        for file in files:
            if '_genomes.' in file and '.fa' in file:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                if str(bin_id) not in selected_bin.keys():
                    n+=1
                    os.system('rm '+file)
    print 'Removed', n, 'refined bins'

    n=0
    os.chdir(pwd+'/'+str(original_bin_folder))
    for root, dirs, files in os.walk(pwd+'/'+str(original_bin_folder)):
        for file in files:     
            if '_genomes.' in file:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                if bin_id in bestbin.keys():
                    if bin_id in selected_bin.keys():
                        n+=1
                        os.system('cp '+file+' '+pwd+'/'+str(new_bin_folder))
                else:
                    n+=1
                    os.system('cp '+file+' '+pwd+'/'+str(new_bin_folder))
                    ids_list=file.split('.')
                    ids_list.remove(ids_list[-1])
                    bin_id='.'.join(ids_list)
                    selected_bin[bin_id]=1
    print 'Moved', n, 'to refined bins folder'

    os.chdir(pwd+'/'+str(new_bin_folder))
    f_bin_checkm=open(new_bin_folder+'_retrieved_bin_stats_ext.tsv', 'w')
    for item in new_bins_checkm.keys():
        if item in selected_bin.keys():
            f_bin_checkm.write(str(item)+'\t'+str(new_bins_checkm[item])+'\n')

    # print selected_bin
    for item in selected_bin.keys():
        if item in bin_checkm.keys():
            # print 'Selected', item
            f_bin_checkm.write(str(item)+'\t'+str(bin_checkm[item])+'\n')
    f_bin_checkm.close()

    os.chdir(pwd)

def parse_bin_in_bestbinset(assemblies_list, binset, BestBinSet_list, PE_connections_list):
    pwd=os.getcwd()
    assemblies={}
    for item in assemblies_list:
        assemblies[item]=[]

    ### Record bins from source
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        os.chdir(pwd+'/'+str(binset))
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('_genomes.')[1]
                if '.fa' in hz:
                    for item in assemblies.keys():
                        if item in file:
                            assemblies[item].append(file)
    os.chdir(pwd)
    os.system('mkdir '+binset+'_recruited')
    kmer_list, bin_connecting_contigs_total, connection_total, bin_kmer_total, kmer_total=[], {}, {}, {}, {}
    for i in range(0, len(assemblies_list)):
        print 'Parsing', str(assemblies_list[i]), 'related bins'
        print assemblies[assemblies_list[i]]
        os.system('perl calc.kmerfreq.pl -i '+str(assemblies_list[i])+' -o '+str(assemblies_list[i])+'.kmer.txt')
        kmer_list.append(str(assemblies_list[i])+'.kmer.txt')
        A=PE_connecting_contigs(assemblies_list[i], PE_connections_list[i], binset, str(assemblies_list[i])+'.kmer.txt')
        bin_connecting_contigs=A[0]
        bin_kmer=A[1]
        contig_kmer=A[2]
        connections_single=A[3]
        bin_connecting_contigs_total.update(bin_connecting_contigs)
        bin_kmer_total.update(bin_kmer)
        kmer_total.update(contig_kmer)
        connection_total.update(connections_single)
        print 'Parsed', str(assemblies_list[i])
        print '--------------------'

    # f=open('Total_kmer.txt','w')
    # for item in kmer_total.keys():
    #     f.write(str(item)+'\t'+str(kmer_total[item])+'\n')
    # f.close()

    # print kmer_total
    f=open('Bin_connecting_contigs.txt', 'w')
    f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
    for item in bin_connecting_contigs_total.keys():
        f.write(str(item)+'\t'+str(bin_connecting_contigs_total[item])+'\n')
    f.close()
    
    A=record_bin_coverage(binset)
    bin_contig_cov=A[0]
    bin_contig=A[1]
    contig_cov=A[2]
    total_bin_contigs=A[3]

    bin_extract_contig, selected_1st, elemimated_contig, elemimated_contig_total, TNFs_exceptional_contigs, m={}, {}, {}, {}, {}, 0
    for bin in bin_connecting_contigs_total.keys():
        m+=1
        n1=0
        print  bin, 'coverage filtration of single copy contigs'
        for connecting_contig in bin_connecting_contigs_total[bin].keys():
            n1+=1
            # print 'Testing', connecting_contig
            test_contig_bin_cov={}
            test_contig_bin_cov.update(bin_contig_cov[bin])
            contig_num=len(bin_contig_cov[bin])

            if connecting_contig in contig_cov.keys():
                test_contig_bin_cov[connecting_contig]=contig_cov[connecting_contig]
                num_coverage=len(contig_cov[connecting_contig])

            ### Filtration of contigs
            cov_index_total_num, cov_index={}, {}
            for item in bin_contig_cov[bin].keys():  
                for i in range(0, num_coverage):
                    if i not in cov_index_total_num.keys():
                        cov_index_total_num[i]=float(bin_contig_cov[bin][item][i+1])
                        cov_index[i]=[]
                        cov_index[i].append(float(bin_contig_cov[bin][item][i+1]))
                    else:
                        cov_index_total_num[i]+=float(bin_contig_cov[bin][item][i+1])
                        cov_index[i].append(float(bin_contig_cov[bin][item][i+1]))
            
            judgement, upper, lower = 0, {}, {}
            for item in cov_index.keys():
                cov_index[item].sort() ### sort low to high
                list_key_num=len(cov_index[item])
                p25=int(0.25*list_key_num)
                p75=int(0.75*list_key_num)

                n=0
                for i in cov_index[item]:
                    n+=1
                    if n == p25:
                        Q1=i
                    elif n == p75:
                        Q3=i
                    else:
                        continue

                IQR = Q3 - Q1
                upper[item+1] = Q3 + 1.5 * IQR
                lower[item+1] = Q1 - 1.5 * IQR
            
            for i in range(1, num_coverage+1):
                if contig_cov[connecting_contig][i] <= upper[i] and contig_cov[connecting_contig][i] >= lower[i]:
                    # print str(i), str(contig_cov[connecting_contig][i]), str(upper[i]), str(lower[i])
                    judgement+=1
                
            if judgement == num_coverage:
                if bin not in selected_1st.keys():
                    selected_1st[bin]={}
                    selected_1st[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    selected_1st[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]

                coverage_data, contigs_ids, coverage_list, n, test_index={}, [], [], 0, 0
                for item in test_contig_bin_cov.keys():
                    n+=1
                    contigs_ids.append(item)
                    for i in range(1, num_coverage+1):
                        if i not in coverage_data.keys():
                            coverage_data[i]=[]
                            # print i
                        coverage_data[i].append(test_contig_bin_cov[item][i])
                        coverage_list.append(test_contig_bin_cov[item][i])
                
                    if item == connecting_contig:
                        test_index=n-1

                coverage_array=np.array(coverage_list).reshape((n,num_coverage))

                A=PCA_slector(coverage_array, n)
                newData=A[0]
                explained_variance_ratio=A[1]
                bin_outlier=test_outlier(connecting_contig, newData, test_index)
                if bin_outlier == 1:
                    if bin not in bin_extract_contig.keys():
                        bin_extract_contig[bin]={}
                        bin_extract_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    else:
                        bin_extract_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    if bin not in elemimated_contig.keys():
                        elemimated_contig[bin]={}
                        elemimated_contig_total[bin]={}
                        elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                        elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    else:
                        elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                        elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            else:
                if bin not in elemimated_contig.keys():
                    elemimated_contig[bin]={}
                    elemimated_contig_total[bin]={}
                    elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]

    bin_extract_contig_TNF, elemimated_contig_TNF={}, {}
    for bin in bin_extract_contig.keys():
        print  bin, 'TNFs filtration of contigs'
        for connecting_contig in bin_extract_contig[bin].keys():
            Bins_TNFs_test=[]
            test_contig_bin_kmer={}
            test_contig_bin_kmer.update(bin_kmer_total[bin])
            test_contig_bin_kmer[connecting_contig]=kmer_total[connecting_contig]
            # print connecting_contig
            # print kmer_total[connecting_contig]
            num_contig=len(test_contig_bin_kmer)
            n1, n_connecting_contig=0, 0
            for item in test_contig_bin_kmer.keys():
                n1+=1
                if item == connecting_contig:
                    n_connecting_contig=n1-1
                lis=test_contig_bin_kmer[item].split('\t')
                for i in range(0, len(lis)):
                    Bins_TNFs_test.append(lis[i])

            # f=open('test.txt','w')
            # for item in Bins_TNFs_test:
            #     f.write(item+'\n')
            # f.close()

            TNF_array=np.array(Bins_TNFs_test).reshape((num_contig, 256))
            
            try:
                # TNF_array.dropna(inplace=True)
                A=PCA_slector(TNF_array, num_contig)
                newData=A[0]
                explained_variance_ratio=A[1]
                bin_outlier=test_outlier(connecting_contig, newData, n_connecting_contig)

                if bin_outlier == 1:
                    if bin not in bin_extract_contig_TNF.keys():
                        bin_extract_contig_TNF[bin]={}
                        bin_extract_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    else:
                        bin_extract_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    if bin not in elemimated_contig_TNF.keys():
                        elemimated_contig_TNF[bin]={}
                        elemimated_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    else:
                        elemimated_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]

                    if bin not in elemimated_contig_total.keys():
                        elemimated_contig_total[bin]={}
                        elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    else:
                        elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            except:
                print bin, 'contig-', connecting_contig, 'TNFs clustering error'
                print 'Adding', connecting_contig, 'to eliminated contigs list'
                if bin not in TNFs_exceptional_contigs.keys():
                    TNFs_exceptional_contigs[bin]={}
                    TNFs_exceptional_contigs[bin][connecting_contig]=kmer_total[connecting_contig]
                else:
                    TNFs_exceptional_contigs[bin][connecting_contig]=kmer_total[connecting_contig]

                if bin not in elemimated_contig_TNF.keys():
                    elemimated_contig_TNF[bin]={}
                    elemimated_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    elemimated_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]

                if bin not in elemimated_contig_total.keys():
                    elemimated_contig_total[bin]={}
                    elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]

    f=open('Coverage_1st_filtrated_bin_connecting_contigs.txt','w')
    for bin in selected_1st.keys():
        f.write(str(bin)+'\t'+str(selected_1st[bin])+'\n')
    f.close()

    f=open('Coverage_2nd_filtrated_bin_connecting_contigs.txt','w')
    for bin in bin_extract_contig.keys():
        f.write(str(bin)+'\t'+str(bin_extract_contig[bin])+'\n')
    f.close()

    f=open('Coverage_eliminated_bin_connecting_contigs.txt','w')
    for bin in elemimated_contig.keys():
        f.write(str(bin)+'\t'+str(elemimated_contig[bin])+'\n')
    f.close()

    f=open('TNF_filtrated_bin_connecting_contigs.txt','w')
    for bin in bin_extract_contig_TNF.keys():
        f.write(str(bin)+'\t'+str(bin_extract_contig_TNF[bin])+'\n')
    f.close()

    f=open('TNF_eliminated_bin_connecting_contigs.txt','w')
    for bin in elemimated_contig_TNF.keys():
        f.write(str(bin)+'\t'+str(elemimated_contig_TNF[bin])+'\n')
    f.close()

    f=open('Total_eliminated_bin_connecting_contigs.txt','w')
    for bin in elemimated_contig_total.keys():
        f.write(str(bin)+'\t'+str(elemimated_contig_total[bin])+'\n')
    f.close()

    f=open('TNFs_exceptional_contigs.txt','w')
    for bins in TNFs_exceptional_contigs.keys():
        for contigs in TNFs_exceptional_contigs[bins].keys():
            f.write(str(bins)+'\t'+str(contigs)+'\t'+str(TNFs_exceptional_contigs[bins][contigs])+'\n')
    f.close()

    print 'Writing retrieved bins'
    os.system('mkdir '+binset+'_retrieved')
    os.system('mv *_filtrated_bin_connecting_contigs.txt *_eliminated_bin_connecting_contigs.txt Bin_connecting_contigs* '+binset+'_retrieved')
    os.chdir(binset+'_retrieved')
    for item in bin_extract_contig.keys():
        bin_name_list=str(item).split('.')
        bin_name_list.remove(bin_name_list[-1])
        bin_name='.'.join(bin_name_list)
        f=open(bin_name+'.retrieved.fa','w')
        if item in bin_contig.keys():
            for ids in bin_contig[item].keys():
                f.write('>'+str(ids)+'\n'+str(bin_contig[item][ids])+'\n')
        
        for contigs in bin_extract_contig[item].keys():
            f.write('>'+str(contigs)+'\n'+str(total_bin_contigs[contigs])+'\n')
        f.close()

    os.chdir(pwd)
    print 'Recalculation of connections file of retrieved bins'
    new_bin_connections={}
    for bin in bin_extract_contig.keys():
        for contigs in bin_extract_contig[bin].keys():
            if contigs in connection_total.keys():
                for connecting_contigs in connection_total[contigs].keys():
                    if connecting_contigs not in  bin_extract_contig[bin].keys() and connecting_contigs not in bin_contig[bin].keys():
                        if bin not in new_bin_connections.keys():
                            new_bin_connections[bin]=int(connection_total[contigs][connecting_contigs])
                        else:
                            new_bin_connections[bin]+=int(connection_total[contigs][connecting_contigs])
        for contigs in bin_contig[bin].keys():
            if contigs in connection_total.keys():
                for connecting_contigs in connection_total[contigs].keys():
                    if connecting_contigs not in  bin_extract_contig[bin].keys() and connecting_contigs not in bin_contig[bin].keys():
                        if bin not in new_bin_connections.keys():
                            new_bin_connections[bin]=int(connection_total[contigs][connecting_contigs])
                        else:
                            new_bin_connections[bin]+=int(connection_total[contigs][connecting_contigs])
                    
    os.chdir(binset+'_retrieved')
    new_bin_connections2, new_bin_connections3={}, {}
    f=open('Retrieved_bins_connections.txt','w')
    for item in new_bin_connections.keys():
        bin_name_list=str(item).split('.')
        bin_name_list.remove(bin_name_list[-1])
        bin_name='.'.join(bin_name_list)
        bin_name_new=bin_name+'.retrieved.fa'
        bin_name_new2=bin_name+'.retrieved'
        f.write(str(bin_name_new)+'\t'+str(new_bin_connections[item])+'\n')
        new_bin_connections2[str(bin_name_new)]=str(new_bin_connections[item])
        new_bin_connections3[str(bin_name_new2)]=str(new_bin_connections[item])
    f.close()
    
    os.chdir(pwd)
    print 'Checking quality of retrieved bins'
    bins_checkm=checkm(str(binset)+'_retrieved')
    for bin in bins_checkm.keys():
        if bin in new_bin_connections3.keys():
            bins_checkm[bin]['Connections']=int(new_bin_connections3[bin])
        else:
            bins_checkm[bin]['Connections']=0

    bin_comparison(str(binset), bins_checkm, str(binset)+'_retrieved')
       
def Contig_recruiter_main(binset):
    pwd=os.getcwd()
    os.chdir(pwd+'/'+str(binset))
    assembly_dict, assembly_list, BestBinSet_list, PE_connections_list={}, [], [], []
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        for file in files:
            if '_genomes.' in file:
                qz=file.split('_genomes.')[0]
                assembly_name_list=qz.split('_')
                assembly_name_list.remove(assembly_name_list[-1])
                assembly_name_list.remove(assembly_name_list[-1])
                assembly_name='_'.join(assembly_name_list)
                assembly_dict[assembly_name]=1
    
    os.chdir(pwd)
    for item in assembly_dict.keys():
        assembly_list.append(item)
        ids_list=item.split('_')
        ids_list.remove(ids_list[0])
        connection_name='_'.join(ids_list)
        PE_connections_list.append('condense_connections_'+connection_name+'.txt')
        BestBinSet_list.append(item+'_BestBinsSet')

    parse_bin_in_bestbinset(assembly_list, binset, BestBinSet_list, PE_connections_list)
    os.chdir(pwd)

if __name__ == '__main__': 
    best_binset_from_multi_assemblies='BestBinset_TNFs_refined'
    Contig_recruiter_main(best_binset_from_multi_assemblies)
