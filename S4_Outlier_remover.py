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

def record_bin_coverage(best_binset_from_multi_assemblies):
    pwd=os.getcwd()
    bin_contigs, bin_contigs_mock, total_bin_contigs, assembly_list, m={}, {}, {}, {}, 0
    print 'Parsing bins'
    for root, dirs, files in os.walk(pwd+'/'+str(best_binset_from_multi_assemblies)):
        os.chdir(pwd+'/'+str(best_binset_from_multi_assemblies))
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
                        total_bin_contigs[str(record.id)]=str(record.seq)
                        bin_contigs[file][str(record.id)]=str(record.seq)
                        bin_contigs_mock[file][str(record.id)]=0
    print 'Parsed', m, 'bins'
    
    os.chdir(pwd)
    f=open('Refined_total_bins_contigs.fa', 'w')
    for item in total_bin_contigs.keys():
        f.write('>'+str(item)+'\n'+str(total_bin_contigs[item])+'\n')
    f.close()

    coverage=open('Coverage_matrix_HBS_refined.txt', 'w')
    coverage.write('Name'+'\t'+'Length'+'\t'+'totalCoverage'+'\t'+'avgCoverage'+'Coverage'+'\n')
    connection_file=open('condense_connections_HBS_refined.txt', 'w')
    connection_file.write('node1'+'\t'+'interaction'+'\t'+'node2'+'\t'+'connections'+'\n')
    for item in assembly_list.keys():
        n=0
        for line in open('Coverage_matrix_for_binning_'+item+'.txt','r'):
            n+=1
            if n >= 2:
                coverage.write(line)

        ids_list=item.split('_')
        ids_list.remove(ids_list[0])
        connection_name='_'.join(ids_list)
        n=0
        for line in open('condense_connections_'+connection_name+'.txt','r'):
            n+=1
            if n >= 2:
                connection_file.write(line)
    coverage.close()
    connection_file.close()

    print 'Recording the coverage of contigs from bins'
    n, contig_cov=0, {}
    for line in open('Coverage_matrix_HBS_refined.txt','r'):
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
    os.system('mkdir bin_coverage')
    os.chdir('bin_coverage')
    for bins in bin_contig_cov.keys():
        f=open(bins+'_coverage_matrix.txt', 'w')
        f.write('Bin'+'\t'+'Coverage'+'\n')
        for contigs in bin_contig_cov[bins].keys():
            f.write(str(contigs)+'\t'+str(bin_contig_cov[bins][contigs])+'\n')
        f.close()
    os.chdir(pwd)
    return bin_contig_cov, bin_contigs, 'Refined_total_bins_contigs.fa', 'condense_connections_HBS_refined.txt'

def outliner_remover(bin_id, contigs_ids, threshold, item_data, explained_variance_ratio, bin_outliner):
    print 'Finding outliner from', bin_id
    four = pd.Series(item_data).describe()
    for item in threshold:
        f1=open('Outlier_in_threshold'+str(item)+'_'+str(bin_id)+'.txt', 'w')
        f2=open('Summary_threshold'+str(item)+'_'+str(bin_id)+'.txt', 'w')
        print(four)
        print('Q1= {0}, Q2= {1}, Q3={2}'.format(four['25%'],four['50%'],four['75%']))
        Q1 = four['25%']
        Q3 = four['75%']
        IQR = Q3 - Q1
        upper1 = Q3 + float(item) * IQR
        lower1 = Q1 - float(item) * IQR
        print(upper1, lower1)
        n1, outliner_record=0, {}
        for i in range(0, len(item_data)):
            if item_data[i] > float(upper1) or item_data[i] < float(lower1):
                n1+=1
                f1.write(str(contigs_ids[i])+'\t'+str(item_data[i])+'\n')
                if str(item) not in bin_outliner[str(bin_id)].keys():
                    bin_outliner[str(bin_id)][str(item)]=[]
                    bin_outliner[str(bin_id)][str(item)].append(str(contigs_ids[i]))
                else:
                    bin_outliner[str(bin_id)][str(item)].append(str(contigs_ids[i]))
        f1.close()
        f2.write(str(four)+'\n'+str('Q1= {0}, Q2= {1}, Q3={2}'.format(four['25%'],four['50%'],four['75%']))+'\n'+'Upper:'+str(upper1)+'\t'+'Lower:'+str(lower1)+'\n'+str(n1)+' outliers in '+str(bin_id)+' under the threshold of '+str(item)+'\n'+'Explained variance ratio:'+str(explained_variance_ratio)+'\n')
        f2.close()
        print n1, 'outliers in', str(bin_id), 'with threshold of', item
        print '-------------------------'
    return bin_outliner

def PCA_slector(data_array, num_contig):
    pca = PCA(n_components=1)
    pca.fit(data_array)
    explained_variance_ratio=pca.explained_variance_ratio_
    print(explained_variance_ratio)
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

def cov_materix(bin_contig_cov, bin_contig, threshold, folder_binset):
    pwd=os.getcwd()
    print 'Transfroming'
    os.system('mkdir bin_cov_outliner')
    bin_outliner={}
    for bin_id in bin_contig_cov.keys():
        bin_outliner[bin_id]={}
        coverage_data, contigs_ids, coverage_list, num_contig={}, [], [], 0
        for contig in bin_contig_cov[bin_id].keys():
            num_contig+=1
            contigs_ids.append(contig)
            num_coverage=len(bin_contig_cov[bin_id][contig])
            for i in range(1, num_coverage+1):
                if i not in coverage_data.keys():
                    coverage_data[i]=[]
                coverage_data[i].append(bin_contig_cov[bin_id][contig][i])
                coverage_list.append(bin_contig_cov[bin_id][contig][i])

        coverage_array=np.array(coverage_list).reshape((num_contig,num_coverage))

        os.chdir('bin_cov_outliner')
        A=PCA_slector(coverage_array, num_contig)
        newData=A[0]
        explained_variance_ratio=A[1]
        bin_outliner=outliner_remover(bin_id, contigs_ids, threshold, newData, explained_variance_ratio, bin_outliner)
        os.chdir(pwd)
    
    os.system('mkdir '+str(folder_binset)+'_coverage_refined')
    os.chdir(str(folder_binset)+'_coverage_refined')
    outliner_sum=open('Summary_ct_outliners.txt','w')
    outliner_sum.write('Bin'+'\t'+'Threshold'+'\t'+'Outliners'+'\n')
    for item in bin_outliner.keys():
        for th in bin_outliner[item].keys():
            if len(bin_outliner[item][th]) != 0:
                outliner_sum.write(str(item)+'\t'+str(th)+'\t'+str(bin_outliner[item][th])+'\n')
    outliner_sum.close()

    for item in bin_outliner.keys():
        dereplicate={}
        for th in bin_outliner[item].keys():
            if len(bin_outliner[item][th]) not in dereplicate.keys():
                dereplicate[len(bin_outliner[item][th])]=1
            else:
                del bin_outliner[item][th]

    for item in bin_contig.keys():
        if item in bin_outliner.keys():
            for th in bin_outliner[item].keys():
                f=open(item+'_ct_'+str(th)+'.fa','w')
                for contig in bin_contig[item].keys():
                    if contig not in bin_outliner[item][th]:
                        f.write('>'+str(contig)+'\n'+str(bin_contig[item][contig])+'\n')
                f.close()
    os.chdir(pwd)

def TNFs_refiner(binset, assembly, coverage_refined_folder, threshold):
    pwd=os.getcwd()
    os.system('mkdir '+binset+'_TNFs_outliner')
    os.system('mkdir '+binset+'_TNFs')
    print 'Calculating TNFs of', assembly
    os.system('perl calc.kmerfreq.pl -i '+str(assembly)+' -o '+str(assembly)+'.kmer.txt')
    
    os.chdir(pwd+'/'+coverage_refined_folder)
    bin_contigs, bin_contigs_mock, bin_outliner={}, {}, {}
    for root, dirs, files in os.walk(pwd+'/'+coverage_refined_folder):
        for file in files:
            if '_genomes.' in file:
                hz=file.split('_genomes.')[1]
                if '.fa' in hz:
                    bin_outliner[file]={}
                    bin_contigs[file]={}
                    bin_contigs_mock[file]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_contigs[file][record.id]=record.seq
                        bin_contigs_mock[file][record.id]=1
    
    os.chdir(pwd+'/'+binset+'_TNFs')
    Bins_TNFs, bin_contig_list={}, {}
    for item in bin_contigs.keys():
        print 'Splitting kmer file,', item
        f=open(item+'_kmer.txt', 'w')
        n, Bins_TNFs[item], bin_contig_list[item] = 0, [], []
        for line in open(pwd+'/'+assembly+'.kmer.txt', 'r'):
            n+=1
            if n >= 2:
                contig=str(line).strip().split('\t')[0]
                lis=str(line).strip().split('\t')
                num=len(bin_contigs_mock[item])
                bin_contigs_mock[item][contig]=1
                if len(bin_contigs_mock[item]) == num:
                # if contig in bin_contigs[item]:
                    f.write(str(line))
                    bin_contig_list[item].append(contig)
                    for i in range(1, len(lis)):
                        Bins_TNFs[item].append(lis[i])
                else:
                    del bin_contigs_mock[item][contig]

    for item in Bins_TNFs.keys():
        print 'Processing TNFs of bin', item
        num_contig=len(bin_contigs[item])
        TNF_array=np.array(Bins_TNFs[item]).reshape((num_contig, 256))
        os.chdir(pwd+'/'+binset+'_TNFs_outliner')
        A=PCA_slector(TNF_array, num_contig)
        newData=A[0]
        explained_variance_ratio=A[1]
        bin_outliner=outliner_remover(item, bin_contig_list[item], threshold, newData, explained_variance_ratio, bin_outliner)
        os.chdir(pwd)
    
    os.system('mkdir '+str(binset)+'_TNFs_refined')
    os.chdir(str(binset)+'_TNFs_refined')
    outliner_sum=open('Summary_TNFs_outliners.txt','w')
    outliner_sum.write('Bin'+'\t'+'Threshold'+'\t'+'Outliners'+'\n')
    for item in bin_outliner.keys():
        for th in bin_outliner[item].keys():
            if len(bin_outliner[item][th]) != 0:
                outliner_sum.write(str(item)+'\t'+str(th)+'\t'+str(bin_outliner[item][th])+'\n')
    outliner_sum.close()

    for item in bin_outliner.keys():
        dereplicate={}
        for th in bin_outliner[item].keys():
            if len(bin_outliner[item][th]) not in dereplicate.keys():
                dereplicate[len(bin_outliner[item][th])]=1
            else:
                del bin_outliner[item][th]

    for item in bin_contigs.keys():
        for th in bin_outliner[item].keys():
            f=open(item+'_TNFs_'+str(th)+'.fa','w')
            for contig in bin_contigs[item].keys():
                if contig not in bin_outliner[item][th]:
                    f.write('>'+str(contig)+'\n'+str(bin_contigs[item][contig])+'\n')
            f.close()
    os.chdir(pwd)

def parse_connections(binset, PE_connections_file):
    print 'Parsing PE-connections file'
    pwd=os.getcwd()
    os.chdir(pwd+'/'+binset)
    bin_connections, bin_contig={}, {}
    for root, dirs, files in os.walk(pwd+'/'+binset):
        for file in files:
            bin_contig[file]=[]
            for record in SeqIO.parse(file, 'fasta'):
                bin_contig[file].append(record.id)

    os.chdir(pwd)
    for item in bin_contig.keys():
        n, connections=0, {}
        for line in open(PE_connections_file,'r'):
            n+=1
            if n >= 2:
                id1=str(line).strip().split('\t')[0]
                id2=str(line).strip().split('\t')[2]
                connections_no=int(str(line).strip().split('\t')[3])
                if id1 in bin_contig[item]:
                    if id2 not in bin_contig[item]:
                        if item not in bin_connections.keys():
                            bin_connections[item]=connections_no
                        else:
                            bin_connections[item]+=connections_no
                elif id2 in bin_contig[item]:
                    if item not in bin_connections.keys():
                        bin_connections[item]=connections_no
                    else:
                        bin_connections[item]+=connections_no
    
    os.chdir(pwd+'/'+binset)
    f=open('Refined_bin_connections.txt', 'w')
    for item in bin_connections.keys():
        f.write(item+'\t'+str(bin_connections[item])+'\n')
    f.close()
    os.chdir(pwd)
    return bin_connections

def bin_comparison(binset, refined_binset, bin_format, bin_connections, com_type):
    pwd=os.getcwd()
    os.system('checkm lineage_wf -t 42 -x fa '+str(refined_binset)+' '+str(refined_binset)+'_checkm')
    os.chdir(refined_binset+'_checkm/storage/')
    refined_checkm={}
    print 'Parsing BestBinset_refined checkm output'
    for line in open('bin_stats_ext.tsv','r'):
        binID=str(line).strip().split('{\'')[0].strip()
        genome_size=str(line).strip().split('Genome size\':')[1].split(', \'Longest')[0].strip()
        taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
        completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
        contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].strip()
        GC=round(float(str(line).strip().split('\'GC\':')[1].split(', \'GCN4\'')[0].strip())*100, 1)
        refined_checkm[str(binID)]={}
        if str(binID)+'.fa' in bin_connections.keys():
            refined_checkm[str(binID)]['Connections']=int(bin_connections[str(binID)+'.fa'])
        else:
            refined_checkm[str(binID)]['Connections']=0
        refined_checkm[str(binID)]['marker lineage']=str(taxon)
        refined_checkm[str(binID)]['Completeness']=float(completeness)
        refined_checkm[str(binID)]['Genome size']=float(genome_size)
        refined_checkm[str(binID)]['Contamination']=float(contamination)

    os.chdir(pwd+'/'+str(binset))
    bin_checkm={}
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        for file in files:
            if '_bin_stats_ext.tsv' in file:
                for line in open(file, 'r'):
                    binID=str(line).strip().split('{\'')[0].strip()
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

    print 'Comparing bins before and after refining process'
    os.chdir(pwd+'/'+str(refined_binset))
    bin_comparison, bin_comparison_list={}, {}
    for refined_bin in refined_checkm.keys():
        if com_type == 'coverage':
            ids=refined_bin.split('_ct_')[0]
        elif com_type == 'TNFs':
            ids=refined_bin.split('_TNFs_')[0]
        ids_list=ids.split('.')
        ids_list.remove(ids_list[-1])
        bin_id='.'.join(ids_list)
        # print refined_bin
        # print bin_id
        if bin_id in bin_checkm.keys():
            if bin_id not in bin_comparison.keys():
                bin_comparison[bin_id]=bin_id+'---'+refined_bin
                bin_comparison_list[bin_id]=[]
                bin_comparison_list[bin_id].append(bin_id)
                bin_comparison_list[bin_id].append(refined_bin)
            else:
                bin_comparison[bin_id]+='---'+refined_bin
                bin_comparison_list[bin_id].append(refined_bin)
    
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
                write_out+='---'+str(refined_checkm[refined_bin_id])
                re_connections=int(refined_checkm[refined_bin_id]['Connections'])
                re_cpn=float(refined_checkm[refined_bin_id]['Completeness'])
                re_ctn=float(refined_checkm[refined_bin_id]['Contamination'])
                re_taxon=str(refined_checkm[refined_bin_id]['marker lineage'])
                re_genome_size=float(refined_checkm[refined_bin_id]['Genome size'])
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
                    bestbin[item][refined_bin_id]=refined_checkm[refined_bin_id]
                elif re_delta == ori_delta:
                    if re_connections < ori_connections:
                        del bestbin[item][ori_bin]
                        bestbin[item][refined_bin_id]=refined_checkm[refined_bin_id]
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
    for root, dirs, files in os.walk(pwd+'/'+str(refined_binset)):
        for file in files:
            if '_genomes.' in file:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                if bin_id not in selected_bin.keys():
                    n+=1
                    os.system('rm '+file)
    print 'Removed', n, 'refined bins'

    n=0
    os.chdir(pwd+'/'+str(binset))
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        for file in files:     
            if '_genomes.' in file:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                if bin_id in bestbin.keys():
                    if bin_id in selected_bin.keys():
                        n+=1
                        os.system('cp '+file+' '+pwd+'/'+str(refined_binset))
                else:
                    n+=1
                    os.system('cp '+file+' '+pwd+'/'+str(refined_binset))
                    ids_list=file.split('.')
                    ids_list.remove(ids_list[-1])
                    bin_id='.'.join(ids_list)
                    selected_bin[bin_id]=1
    print 'Moved', n, 'to refined bins folder'

    os.chdir(pwd+'/'+refined_binset)
    f_bin_checkm=open(refined_binset+'_bin_stats_ext.tsv', 'w')
    for item in refined_checkm.keys():
        if item in selected_bin.keys():
            f_bin_checkm.write(str(item)+'\t'+str(refined_checkm[item])+'\n')

    for item in bin_checkm.keys():
        if item in selected_bin.keys():
            f_bin_checkm.write(str(item)+'\t'+str(bin_checkm[item])+'\n')
    f_bin_checkm.close()

    os.chdir(pwd)

def contig_outlier_remover_main(best_binset_from_multi_assemblies):
    threshold=['1.5', '2', '2.5', '3']

    ### Coverage refine process
    A=record_bin_coverage(best_binset_from_multi_assemblies)
    bin_contig_cov=A[0]
    bin_contig=A[1]
    assembly=A[2]
    PE_connections_file=A[3]
    cov_materix(bin_contig_cov, bin_contig, threshold, best_binset_from_multi_assemblies)
    bin_connections=parse_connections(best_binset_from_multi_assemblies+'_coverage_refined', PE_connections_file)
    bin_comparison(best_binset_from_multi_assemblies, best_binset_from_multi_assemblies+'_coverage_refined', 'fa', bin_connections, 'coverage')
    ### TNFs refine process
    TNFs_refiner(best_binset_from_multi_assemblies, assembly, best_binset_from_multi_assemblies+'_coverage_refined', threshold)
    bin_connections=parse_connections(best_binset_from_multi_assemblies+'_TNFs_refined', PE_connections_file)
    bin_comparison(best_binset_from_multi_assemblies+'_coverage_refined', best_binset_from_multi_assemblies+'_TNFs_refined', 'fa', bin_connections, 'TNFs')

if __name__ == '__main__': 
    best_binset_from_multi_assemblies='BestBinset'
    contig_outlier_remover_main(best_binset_from_multi_assemblies)
