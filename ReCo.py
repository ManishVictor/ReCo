import os				#importing all the required packages
from Bio import SeqIO
import random
import subprocess
import shutil
class corec():
    def cr(self,path):
        print('ReCo.py\ndoi:')
        opath=('/'.join(path.split('/')[:-1]))+'/'+'CoRec'
        os.mkdir(opath)
        os.mkdir(os.path.join(opath,'Sequences'))
        pwd=os.getcwd()
        try:
            shutil.rmtree(os.path.join(pwd,'output','CORETREE'))
        except OSError:
            pass
        os.mkdir(os.path.join(opath,'CORETREE'))
        files=os.listdir(path)
        for each in files:
            cmd=('java -jar UBCG.jar extract -bcg_dir '+os.path.join(opath,'Sequences')+' -i '+os.path.join(path,each)+' -label '+each.split('.')[0])+'>'+os.path.join(opath,'log.txt')
            subprocess.call(cmd,shell=True)
            with open(os.path.join(opath,'log.txt'),'r') as file1:
                rdf=file1.readlines()
            for every in rdf:
                with open(os.path.join(opath,'LOG.txt'),'a') as file2:
                    file2.write(every)
        os.remove(os.path.join(opath,'log.txt'))
        cmnd=('java -jar UBCG.jar align -bcg_dir '+os.path.join(opath,'Sequences')+' -prefix CORETREE > '+os.path.join(opath,'align.txt'))
        subprocess.call(cmnd,shell=True)
        with open(os.path.join(opath,'align.txt'),'r') as file3:
            redf=file3.readlines()
        for ech in redf:
                with open(os.path.join(opath,'ALIGN_LOG.txt'),'a') as file4:
                    file4.write(ech)
        os.remove(os.path.join(opath,'align.txt'))
        with open(os.path.join(opath,'ALIGN_LOG.txt'),'r') as file5:
            rdl=file5.readlines()
        ll=(rdl[-2].split(' ')[-1].rstrip('\n')).split("'")[1]
        shutil.copy(os.path.join(pwd,ll),os.path.join(opath,'CORETREE','core.nwk'))
        return(os.path.join(opath,'Sequences'),opath)
    def seqEx(self,patha,pathb):
        print('ReCo.py\ndoi:')
        os.mkdir(os.path.join(pathb,'SeqEx'))
        listoffiles=os.listdir(patha)
        print('Extracting the nucleotide and protein sequences.')
        for t in range(len(listoffiles)):
            with open(os.path.join(patha,listoffiles[t]),'r') as file1:
                rd_file1=file1.readlines()
            new_dat=(rd_file1[0].split('{"data":{"')[1].split(']],"'))
            for i in range(len(new_dat)):
                with open(os.path.join(pathb,'SeqEx',listoffiles[t].split('.bcg')[0]),'a') as file2:
                    file2.write(str(new_dat[i])+'\n')
        print('Number of sequences found',t+1)
        return(pathb,(os.path.join(pathb,'SeqEx')))
    def seqA(self,patha,pathb):
        print('ReCo.py\ndoi:')
        lof=os.listdir(pathb)
        os.mkdir(os.path.join(patha,'seqA'))
        print('Filtering the sequences on the basis of zero found and Paralogs.')
        for u in range(len(lof)):
            with open(os.path.join(pathb,lof[u]),'r') as file3:
                rd_file3=file3.readlines()
            for each in rd_file3:
                count=0
                spech=each.split(':')[1].split(',')[0]
                if(spech == '[0]'):
                    tobwrtn=(each.split(','))
                    for eech in tobwrtn:
                        if('[0]' in eech):
                            count+=1
                            wrtn=','.join(each.split(',')[count:])
                            with open (os.path.join(patha,'seqA',lof[u]),'a') as file4:
                                file4.write('>'+wrtn)
                else:
                    with open (os.path.join(patha,'seqA',lof[u]),'a') as file4:
                        file4.write('>'+each)
        print('Number of sequences found',u+1)
        return(os.path.join(patha,'seqA'))
    def protf(self,patha,pathb):
        print('ReCo.py\ndoi:')
        listof=os.listdir(pathb)
        os.mkdir(os.path.join(patha,'prot'))
        print('Making the final Protein Sequences.')
        for v in range(len(listof)):
            with open(os.path.join(pathb,listof[v]),'r') as file5:
                rd_file5=file5.readlines()
            for every in rd_file5:
                many=[]
                manyseq=[]
                speve=every.split(':')[1].split(',')[0]
                if(speve=='[1'):
                    title=every.split(':')[0].strip('"')
                    seq=every.split(':')[1].split(',')[3].strip('"')
                else:
                    try:
                        speve_sp=(speve.split('[')[1])
                        for w in range(int(speve_sp)-1):
                            if(w==0):
                                many.append(every.split(':')[0].strip('"'))
                            else:
                                many.append(every.split(':')[1].split(',')[2*w+1].strip('"'))
                        title=many[0]
                        for x in range(1,len(many)):
                            if(many[x].startswith('A')):
                                manyseq.append(many[x])
                            try:
                                seq=random.choice(manyseq)
                            except IndexError:
                                continue
                    except ValueError:
                        continue
                with open(os.path.join(patha,'prot',listof[v]+'.fasta'),'a') as file6:
                    file6.write(title+'\n'+seq+'\n')
        return(os.path.join(patha,'prot'))
    def cmnco(self,patha,pathb):
        print('ReCo.py\ndoi:')
        listof=os.listdir(patha)
        seqid=[]
        finalid=[]
        for lineach1 in listof:
            for lineachseq in SeqIO.parse(os.path.join(patha,lineach1),'fasta'):
                idee=(lineachseq.id)
                seqid.append(idee)
        setseq=set(seqid)
        for setid in setseq:
            if(seqid.count(setid) == len(listof)):
                finalid.append(setid)
        ############  Making sequences from the common set of core genes   #############################################
        os.mkdir(os.path.join(pathb,'common_core'))
        for org in listof:
            for orgseq in SeqIO.parse(os.path.join(patha,org),'fasta'):
                for lnum in range(len(finalid)):
                    if(finalid[lnum] in orgseq.id):
                        with open(os.path.join(pathb,'common_core',org),'a') as fileseq:
                            fileseq.write('>'+orgseq.id+'\n'+str(orgseq.seq)+'\n')
        return(os.path.join(pathb,'common_core'))
    def konkat(self,patha,pathb):
        print('ReCo.py\ndoi:')
        print('Concatenating sequences')
        os.mkdir(os.path.join(patha,'for_phipack'))
        filess=os.listdir(pathb)
        for d in range(len(filess)):
            cont=0
            with open(os.path.join(pathb,filess[d]),'r') as file2:
                rd_f2=file2.readlines()
            for q in range(0,len(rd_f2),2):
                with open(os.path.join(patha,'for_phipack','rec'+str(cont)+'.fasta'),'a') as file3:
                    file3.write(rd_f2[q]+rd_f2[q+1])
                cont+=1
        return(os.path.join(patha,'for_phipack'))
    def msal(self,patha,pathb):
        print('ReCo.py\ndoi:')
        print('Performing MSA, Pre-step for finding the recombinant genes')
        os.mkdir(os.path.join(patha,'phipack_mafft'))
        recolist=os.listdir(pathb)
        for evch in recolist:
            align='mafft --quiet '+os.path.join(pathb,evch)+' > '+os.path.join(patha,'phipack_mafft',evch.split('.fasta')[0]+'.aln.fasta')
            os.system(align)
        return(os.path.join(patha,'phipack_mafft'))
    def recnal(self,patha,pathb):
        print('ReCo.py\ndoi:')
        print('Recombination Analysis')
        os.mkdir(os.path.join(patha,'phipack'))
        philist=os.listdir(pathb)
        for phieach in philist:
            phicmd='/usr/bin/phipack-phi -f '+os.path.join(patha,'phipack_mafft',phieach)+' -t A -w 10 -v -g i > '+os.path.join(patha,'phipack',phieach.split('.aln.fasta')[0])
            try:
                os.system(phicmd)
            except ERROR:
                continue
        ################################ extracting recombinant data ##################################
        print('Creating List of Recombinant and Non-Recombinant genes: Tables,Fasta files\n\n')
        val_nrc=[]
        Val_nrc={}
        names_gene=[]
        os.mkdir(os.path.join(patha,'Recombinant-data'))
        nrclist=os.listdir(os.path.join(patha,'phipack'))
        for nreach in range(len(nrclist)):
            with open(os.path.join(patha,'phipack',nrclist[nreach]),'r') as nrcf:
                nrcfl=nrcf.readlines()
            for inr in range(len(nrcfl)):
                if('PHI (Normal):' in nrcfl[inr]):
                    elnrcfl=nrcfl[inr]
                    val_nrc.append(nrclist[nreach]+':'+(elnrcfl.rstrip('\n').split(':'))[1].strip(' '))
        for rnum in range(len(val_nrc)):
            nvar='rec'+str(rnum)
            for start in val_nrc:
                if(nvar == start.split(':')[0]):
                    try:
                        Val_nrc.update({start.split(':')[0]:float(start.split(':')[1])})
                    except ValueError:
                        Val_nrc.update({start.split(':')[0]:'--'})
                    with open(os.path.join(patha,'Recombinant-data','sorted_List_Phi_Values.txt'),'a') as wrfn:
                        wrfn.write(start+'\n')
                        break
        for key, value in Val_nrc.items():
            if(value != '--' and float(value)<0.05):
                names_gene.append(int(key[3:]))
                with open(os.path.join(patha,'Recombinant-data','all_core_genes_having_recombination.txt'),'a') as rrfn:
                    rrfn.write(key+':'+str(value)+'\n')
        ################################ creating set of non-recombinant core genes #####
        serial=[]
        with open(os.path.join(patha,'Recombinant-data','all_core_genes_having_recombination.txt'),'r') as rcmbnf:
            r_rcmbnf=rcmbnf.readlines()
        for evech in r_rcmbnf:
            vaal=evech.split('\n')[0].split(':')[0].split('rec')[1]
            serial.append(vaal)
        return(os.path.join(patha,'Recombinant-data','all_core_genes_having_recombination.txt'))
    def recm(self,patha,pathb,pathc,pathd):
        print('ReCo.py\ndoi:')
        serial=[]
        heeder=[]
        with open(patha,'r') as rcmbnf:
            r_rcmbnf=rcmbnf.readlines()
        for evech in r_rcmbnf:
            vaal=evech.split('\n')[0].split(':')[0].split('rec')[1]
            serial.append(vaal)
        #print(serial)
        for echnum in serial:
            rekseq='rec'+echnum+'.fasta'
            print(rekseq)
            with open(os.path.join(pathb,rekseq),'r') as rekf:
                rd_rekf=rekf.readlines()[0]
            heeder.append(rd_rekf.rstrip('\n'))   
        ########################## making non-recombinant core genes ######################
        os.mkdir(os.path.join(pathc,'for_MSA'))
        listof=os.listdir(pathd)
        for ccgech in listof:
            seqmsa=''
            for ccgseq in SeqIO.parse(os.path.join(pathd,ccgech),'fasta'):
                if(ccgseq.id not in heeder):
                    seqmsa+=ccgseq.seq
            with open(os.path.join(pathc,'for_MSA','all_for_msa.fasta'),'a') as msaf:
                msaf.write('>'+ccgech.split('.fasta')[0]+'\n'+str(seqmsa)+'\n')
        return(os.path.join(pathc,'for_MSA','all_for_msa.fasta'),os.path.join(pathc,'for_MSA'))
############ Main ############################## Main ######################################## Main ##########################################
if __name__=='__main__':
    obj=corec()
    direc=input('Directory name/Directory path of the Scaffold/Contigs')
    (o1path,oo1path)=obj.cr(direc)# pathof scaffold
    (patha,pathb)=obj.seqEx(o1path,oo1path)
    o2path=obj.seqA(patha,pathb)
    o3path=obj.protf(patha,o2path)
    o4path=obj.cmnco(o3path,patha)
    o5path=obj.konkat(patha,o4path)
    o6path=obj.msal(patha,o5path)
    o7path=obj.recnal(patha,o6path)
    (o8path1,o8path2)=obj.recm(o7path,o5path,patha,o4path)
    print('MSA for Non-Recombinant genes')
    os.mkdir(os.path.join(patha,'newick'))
    os.system('mafft --quiet '+o8path1+' > '+os.path.join(o8path2,'all_for_msa.aln.fasta'))
    print('Making Phylogenetic Tree for Non-Recombinant genes')
    os.system('fasttree < '+os.path.join(o8path2,'all_for_msa.aln.fasta')+' > '+os.path.join(patha,'newick','nreco.nwk'))
    print('Results saved at:',oo1path)
