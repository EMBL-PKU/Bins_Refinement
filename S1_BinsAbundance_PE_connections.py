#!/usr/bin/python
# -*- coding: UTF-8 -*-
#coding=utf-8

from Bio import SeqIO
import sys, os
from collections import Counter

def intervalue(Xmin, Xmax, Y, Z):
    delta=Xmax-Xmin
    tim=int(delta/Y)
    for i in range(1, tim+2):
        if Z > Xmin:
            Xmin += Y
        else:
            return Xmin

def covrange(X):
    if float(X)==0:
        return '0000'
    elif float(X) > 0 and  float(X) < 9:
        X+=1
        return '000'+str(int(X))
    elif float(X) == 9:
        return '0009'
    elif float(X) > 9 and  float(X) < 10:
        return '0010'
    elif float(X) >= 10 and float(X) < 20:
        return '00'+str(intervalue(10, 20, 2, X))
    elif float(X) >= 20 and float(X) < 50:
        return '00'+str(intervalue(20, 50, 3, X))
    elif float(X) >= 50 and float(X) < 100:
        return '00'+str(intervalue(50, 100, 5, X))
    elif float(X) >= 100 and  float(X) < 300:
        return '0'+str(intervalue(100, 300, 10, X))
    elif float(X) >= 300 and  float(X) < 700:
        return '0'+str(intervalue(300, 700, 20, X))
    elif float(X) >= 700 and  float(X) < 1000:
        return '0'+str(intervalue(700, 1000, 30, X))
    elif float(X) >= 1000 and  float(X) < 2000:
        return str(intervalue(1000, 2000, 50, X))
    elif float(X) >= 2000 and  float(X) < 5000:
        return str(intervalue(2000, 5000, 100, X))
    else:
        return '10000'

def dcovrange(X):
    if float(X)==0:
        return '0000'
    elif float(X) > 0 and float(X) < 9:
        X+=1
        return '000'+str(int(X))
    elif float(X) == 9:
        return '0009'
    elif float(X) > 9 and float(X) < 20:
        X+=1
        return '00'+str(int(X))
    elif float(X) >= 20 and float(X) < 100:
        return '00'+str(intervalue(20, 100, 2, X))
    elif float(X) >= 100 and  float(X) < 200:
        return '0'+str(intervalue(100, 200, 2, X))
    elif float(X) >= 200 and  float(X) < 1000:
        return '0'+str(intervalue(200, 1000, 5, X))
    elif float(X) >= 1000 and  float(X) < 2000:
        return str(intervalue(1000, 2000, 10, X))
    elif float(X) >= 2000 and  float(X) < 5000:
        return str(intervalue(2000, 5000, 50, X))
    else:
        return '10000'

def CoverageMatrix(depth_file, assembly_name):
    path=os.getcwd()

    fout=open('Coverage_matrix_for_binning_'+str(assembly_name)+'.txt', 'w')

    n, title, cov=0, {}, {}
    for line in open(str(depth_file), 'r'):
        n+=1
        if n == 1:
		    num_cov_groups=int(str(line).strip().count(".bam-var"))+2
		    title['Name']='Length'+'\t'+'totalCoverage'+'\t'+'avgCoverage'
		    for i in range(2, num_cov_groups):
			    m=title['Name']
			    title['Name']=m+'\t'+'Coverage'+str(i-1)+'\t'+'Cov'+str(i-1)+'range'+'\t'+'Cov'+str(i-1)+'drange'
			    if i == 2:
				    covgs='Length'+'\t'+'totalCoverage'+'\t'+'avgCoverage'+'\t'+'Coverage1'+'\t'+'Cov1range'+'\t'+'Cov1drange'
			    else:
				    covgs+='\t'+'Coverage'+str(i-1)+'\t'+'Cov'+str(i-1)+'range'+'\t'+'Cov'+str(i-1)+'drange'
		    fout.write('Name'+'\t'+str(title['Name'])+'\n')
        else:
            totalcov=str(line).strip().split('\t')[2]
            num=float(num_cov_groups)-2
            avgcov=round(float(totalcov)/float(num), 3)
            cov[str(line).strip().split('\t')[0]]=str(line).strip().split('\t')[1]+'\t'+str(totalcov)+'\t'+str(avgcov)
            for i in range(2, num_cov_groups):
                m=cov[str(line).strip().split('\t')[0]]
                covi=str(line).strip().split('\t')[i*2-1]
                covr=covrange(float(covi))
                covdr=dcovrange(float(covi))
                cov[str(line).strip().split('\t')[0]]=m+'\t'+str(covi)+'\t'+str(covr)+'\t'+str(covdr)

    for item in cov.keys():
	    fout.write(str(item)+'\t'+str(cov[str(item)])+'\n')
    fout.close()
    return cov, covgs, 'Coverage_matrix_for_binning_'+str(assembly_name)+'.txt'

def BinAbundance(depth, cov, covgs, output_format, name_of_the_binning_project):
    path=os.getcwd()
    os.chdir(path+'/'+name_of_the_binning_project+'_genomes')
    genome_contig, genome_contig2={}, {}
    for root,dirs,files in os.walk(path+'/'+name_of_the_binning_project+'_genomes'):
        for file in files:
            # print 'Reading', file
            hz=str(file).split('.')[-1]
            if str(output_format) == hz and str(depth) not in str(file):
                file_name_list=str(file).split('.')
                file_name_list.remove(file_name_list[-1])
                file_name='.'.join(file_name_list)
                fout=open(str(file_name)+'_contigs_summary.txt', 'w')
                fout.write('Name'+'\t'+str(covgs)+'\t'+'GC%'+'\n')
                for record in SeqIO.parse(str(file),'fasta'):
                    gc=int(str(record.seq).count("G"))+int(str(record.seq).count("C"))
                    gc_ratio=round(float(100*gc/float(len(record.seq))), 1)
                    fout.write(str(record.id)+'\t'+str(cov[str(record.id)])+'\t'+str(gc_ratio)+'%'+'\n')
                    genome_contig[str(record.id)]=str(file_name)
                    genome_contig2[str(record.id)]=str(file_name)
                fout.close()
            else:
                continue
            
    genome_summary={}
    print '----------------------------------------------------'
    print 'Checking ', name_of_the_binning_project, 'summary file'
    for root,dirs,files in os.walk(path+'/'+name_of_the_binning_project+'_genomes'):
        for file in files:
            if '_contigs_summary.txt' in str(file):
                total_cov_bin=0
                genomeID=str(file).split('_contigs_summary.txt')[0]
                # print 'Checking', file
                n=0
                for line in open(str(file), 'r'):
                    n+=1
                    if n>=2:
                        total_cov_bin+=float(str(line).strip().split('\t')[3])
                    else:
                        continue
                if n>= 2:
                    genome_summary[str(genomeID)]=float(total_cov_bin)/float(n-1)
                #else:
                    #genome_summary[str(genomeID)]=0
            else:
                continue

    gen_sum=sorted(genome_summary.items(), key=lambda genome_summary:genome_summary[1])

    fout=open('Genome_summary_'+str(name_of_the_binning_project)+'.txt', 'w')
    fout.write('Prebin'+'\t'+'PreviousID'+'\t'+'avgCov'+'\n')
    genome_id, genome_avgcov, n={}, {}, 0
    for item in gen_sum:
        n+=1
        genome_id[str(item).split('\'')[1]]=n
        genome_avgcov[str(item).split('\'')[1]]=str(item).split(',')[1].split(')')[0].strip()
        fout.write(str(n)+'\t'+str(item).split('\'')[1].strip()+'\t'+str(item).split(',')[1].split(')')[0].strip()+'\n')
    fout.close()

    for item in genome_contig.keys():
        if  genome_contig[str(item)] in genome_id.keys():
            m=str(genome_contig[str(item)])
            genome_contig[str(item)]=str(genome_id[str(genome_contig[str(item)])])+'\t'+str(m)+'\t'+str(genome_avgcov[str(genome_contig[str(item)])])+'\t'+'---'

    fout=open('Genome_contig_summary_'+str(name_of_the_binning_project)+'.txt', 'w')
    fout.write('ID'+'\t'+'Prebinid'+'\t'+'PreviousID'+'\t'+'avgCov'+'\t'+'EssCompleteness'+'\n')
    for item in genome_contig.keys():
        fout.write(str(item)+'\t'+str(genome_contig[str(item)])+'\n')
    fout.close()

    os.chdir(path+'/'+name_of_the_binning_project+'_genomes')
    # os.system('pwd')
    f=open('Bins_change_ID_'+name_of_the_binning_project+'.txt', 'w')
    try:
        os.mkdir('Original_bins')
    except:
        print 'Folder Original_bins already exist'

    for item in genome_id.keys():
        if 'metabat' in item:
            os.system('mv '+item+'.fa Original_bins/')
            f.write(item+'.fa changed to '+name_of_the_binning_project+'_genomes.'+str(genome_id[item])+'.fa'+'\n')
        elif 'maxbin2' in item or 'concoct' in item:
            os.system('mv '+item+'.fasta Original_bins/')
            f.write(item+'.fa changed to '+name_of_the_binning_project+'_genomes.'+str(genome_id[item])+'.fasta'+'\n')
        else:
            continue
    f.close()
    
    for item in genome_id.keys():
        if 'metabat' in item:
            os.system('cp Original_bins/'+item+'.fa '+name_of_the_binning_project+'_genomes.'+str(genome_id[item])+'.fa')
        elif 'maxbin2' in item or 'concoct' in item:
            os.system('cp Original_bins/'+item+'.fasta '+name_of_the_binning_project+'_genomes.'+str(genome_id[item])+'.fasta')
        else:
            continue

    os.system('mv *_contigs_summary.txt Original_bins/')
    os.system('tar zcvf Original_bins.tar.gz Original_bins')
    os.system('rm -rf Original_bins')

    checkm={}
    print 'Reading', name_of_the_binning_project, 'checkm output'
    f_bin_checkm=open(name_of_the_binning_project+'_bin_stats_ext.tsv', 'w')
    for line in open(str(path)+'/'+str(name_of_the_binning_project)+'_checkm/storage/bin_stats_ext.tsv','r'):
        binID=str(line).strip().split('{\'')[0].strip()
        genome_size=str(line).strip().split('Genome size\':')[1].split(', \'Longest')[0].strip()
        taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
        completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
        contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].strip()
        GC=round(float(str(line).strip().split('\'GC\':')[1].split(', \'GCN4\'')[0].strip())*100, 1)
        checkm[str(binID)]=str(GC)+'\t'+str(genome_size)+'\t'+str(taxon)+'\t'+str(completeness)+'\t'+str(contamination)       
        if str(binID) in genome_id.keys():
            line=line.replace(str(binID), name_of_the_binning_project+'_genomes.'+str(genome_id[binID]))
            f_bin_checkm.write(line)

    for item in gen_sum:
        if str(item).split('\'')[1] not in checkm.keys():
            checkm[str(item).split('\'')[1]]='0'+'\t'+'0'+'\t'+'root'+'\t'+'0'+'\t'+'0'
            f_bin_checkm.write(name_of_the_binning_project+'_genomes.'+str(genome_id[str(item).split('\'')[1]])+'\t'+'{\'GC std\': 0, \'GCN4\', \'Genome size\': 0, \'Longest, \'marker lineage\': \'root\', \'Completeness\': 0.0, \'Contamination\': 0.0}'+'\n')
    f_bin_checkm.close()
    os.system('cp '+name_of_the_binning_project+'_bin_stats_ext.tsv '+str(path)+'/'+str(name_of_the_binning_project)+'_checkm/storage/')

    # os.mkdir('Completeness_lower_10percent_bins')

    fout=open('prebinned_genomes_output_for_dataframe_'+str(name_of_the_binning_project)+'.txt', 'w')
    fout.write('ID'+'\t'+'Prebinid'+'\t'+'PreviousID'+'\t'+'avgCov'+'\t'+'EssCompleteness'+'\t'+'AvgGC'+'\t'+'GenomeSize'+'\t'+'Taxon'+'\t'+'Completeness'+'\t'+'Contamination'+'\n')
    for item in genome_contig.keys():
        fout.write(str(item)+'\t'+name_of_the_binning_project+'_genomes.'+str(genome_contig[str(item)])+'\t'+str(checkm[str(genome_contig2[str(item)])])+'\n')

    for item in cov.keys():
        if item not in genome_contig.keys():
            fout.write(str(item)+'\t'+name_of_the_binning_project+'_genomes.0'+'\t'+'unclustered'+'\t'+'---'+'\t'+'---'+'\t'+'---'+'\t'+'---'+'\t'+'---'+'\t'+'uc'+'\t'+'uc'+'\n')
    fout.close()
    os.chdir(path)
    return 'prebinned_genomes_output_for_dataframe_'+str(name_of_the_binning_project)+'.txt'

def GenerationOfGenomeGroupList(prebin_dataframe, PE_connection_file, name_of_the_binning_project):
    print '---------------------------'
    print 'Reading PE-connections file'

    pwd=os.getcwd()
    os.chdir(pwd+'/'+name_of_the_binning_project+'_genomes')

    contig_genome, m={}, 0
    for line in open(str(prebin_dataframe), 'r'):
        m+=1
        if m >= 2:
            contig_genome[str(line).strip().split('\t')[0]]=str(line).strip().split('\t')[1]

    genome_connection, m={}, 0
    for line in open(pwd+'/'+str(PE_connection_file), 'r'):
        m+=1
        if m >= 2:
            node1=str(line).strip().split('\t')[0]
            node2=str(line).strip().split('\t')[2]
            num_connections=str(line).strip().split('\t')[3]
            # print str(node1), str(node2)
            if str(node1) in contig_genome.keys() and str(node2) in contig_genome.keys() and contig_genome[str(node1)] != contig_genome[str(node2)]:
                if contig_genome[str(node1)] not in genome_connection.keys():
                    genome_connection[contig_genome[str(node1)]]={}
                    genome_connection[contig_genome[str(node1)]][contig_genome[str(node2)]]=str(num_connections)
                else:
                    if contig_genome[str(node2)] not in genome_connection[contig_genome[str(node1)]].keys():
                        genome_connection[contig_genome[str(node1)]][str(contig_genome[str(node2)])]=str(num_connections)
                    else:
                        m1=genome_connection[contig_genome[str(node1)]][str(contig_genome[str(node2)])]
                        genome_connection[contig_genome[str(node1)]][contig_genome[str(node2)]]=str(int(m1)+int(num_connections))
        
                if contig_genome[str(node2)] not in genome_connection.keys():
                    genome_connection[contig_genome[str(node2)]]={}
                    genome_connection[contig_genome[str(node2)]][contig_genome[str(node1)]]=str(num_connections)
                else:
                    if contig_genome[str(node1)] not in genome_connection[contig_genome[str(node2)]].keys():
                        genome_connection[contig_genome[str(node2)]][str(contig_genome[str(node1)])]=str(num_connections)
                    else:
                        m1=genome_connection[contig_genome[str(node2)]][str(contig_genome[str(node1)])]
                        genome_connection[contig_genome[str(node2)]][str(contig_genome[str(node1)])]=str(int(m1)+int(num_connections))

    # print str(genome_connection)

    genome_group='Genome_group_for_cytoscape_'+str(name_of_the_binning_project)+'.txt'
    genome_group_list='Genome_group_all_list_'+str(name_of_the_binning_project)+'.txt'
    fout=open(str(genome_group), 'w')
    fout2=open(str(genome_group_list), 'w')
    # fout3=open('Genome_group_all_top'+str(topGenome)+'.txt', 'w')
    fout.write('Genome1'+'\t'+'Connections'+'\t'+'Genome2'+'\n')
    fout2.write('Genome'+'\t'+'Connecting genomes'+'\n')
    for item in genome_connection.keys():
        fout2.write(str(item)+'\t'+str(genome_connection[item]).strip()+'\n')
        num=len(genome_connection[item])
        if num == 1:
            fout.write(str(item)+'\t'+str(genome_connection[item]).split(':')[1].split('\'')[1].strip()+'\t'+str(genome_connection[item]).split(':')[0].split('\'')[1].strip()+'\n')        
        else:
            lis=str(genome_connection[item]).split(',')
            for i in range(0, num):
                fout.write(str(item)+'\t'+str(lis[i]).split(':')[1].split('\'')[1].strip()+'\t'+str(lis[i]).split(':')[0].split('\'')[1].strip()+'\n')

    fout.close()
    fout2.close()

    genome_total_connection={}
    genome_total_connection_file='Bins_total_connections_'+str(name_of_the_binning_project)+'.txt'
    f=open(str(genome_total_connection_file), 'w')
    f.write('Bin'+'\t'+'Total_connections'+'\n')
    for item in genome_connection.keys():
        genome_total_connection[item]=0
        if len(genome_connection[item]) != 0:
            for i in genome_connection[item].keys():
                genome_total_connection[item]+=int(genome_connection[item][i])
        f.write(str(item)+'\t'+str(genome_total_connection[item])+'\n')
    f.close()

    print 'Generation of Genome Group of',str(name_of_the_binning_project),'List Done!'
    os.chdir(pwd)

    return genome_total_connection_file

def binsabundance_pe_connections(bins_folders_name_list, depth_file, PE_connections_file, assembly_name, coverage_matrix_list):
    coverage_matrix=CoverageMatrix(depth_file, assembly_name)
    coverage_matrix_list.append(coverage_matrix[2])

    for item in bins_folders_name_list:
        if 'maxbin2' in item or 'concoct' in item:
            a=BinAbundance(depth_file, coverage_matrix[0], coverage_matrix[1], 'fasta', item)
            GenerationOfGenomeGroupList(a, PE_connections_file, item)
        elif 'metabat' in item:
            a=BinAbundance(depth_file, coverage_matrix[0], coverage_matrix[1], 'fa', item)
            GenerationOfGenomeGroupList(a, PE_connections_file, item)
        else:
            continue
    return coverage_matrix_list
            
if __name__ == '__main__': 


    bins_folders_name_list=['3_ZH_MO_final.contigs.fa_0.3_maxbin2', '3_ZH_MO_final.contigs.fa_0.5_maxbin2', '3_ZH_MO_final.contigs.fa_0.7', '3_ZH_MO_final.contigs.fa_0.9_maxbin2', '3_ZH_MO_final.contigs.fa_200_metabat', '3_ZH_MO_final.contigs.fa_300_metabat', '3_ZH_MO_final.contigs.fa_400_metabat', '3_ZH_MO_final.contigs.fa_500_metabat', '3_ZH_MO_final.contigs.fa_200_concoct', '3_ZH_MO_final.contigs.fa_400_concoct']
    depth_file='3_assembly.depth.txt'
    PE_connections_file='condense_connections_ZH_MO_final.contigs.fa.txt'
    assembly_name='3_ZH_MO_final.contigs.fa'

    coverage_matrix_list=[]
    binsabundance_pe_connections(bins_folders_name_list, depth_file, PE_connections_file, assembly_name, coverage_matrix_list)
