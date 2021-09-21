# ReCo
Recombination in House-Keeping Genes
# ReCo.py
Recombination Analysis (Linux-like systems)

**Introduction**
House-Keeping gens are the conserved genes or the constitutive genes which are required for the basal cellular function.
Under any condition they are expressed for the maintenance of the cellular functions.Eventhough they show high level of 
conservation it is intuitive to think that they would also show changes over time as evolution is ubiquitous. Genetic 
variability occurs either through mutations or genetic recombinations. For the bacterial house-keeping genes this program 
finds the recombinantion. 
Generally, recombination is related to horizontal gene transfer in bacteria hence this program is helpful for users who 
are studying the evolution of bacteria with respect to the changes in the house-keeping genes.(Gain or loss of funtion).It 
also produces the core genome based phylogenetic tree and core non-recombinant genome based phylogenetic tree.
This program specifically focuses on the housekeeping gens and thus all its results are pertaining to highly conserved genes.

############# This program requires SCAFFOLDS/CONTIGS of the sequences not the annotated sequences  ####################

**How to use it?**
1. Install UBCG (https://www.ezbiocloud.net/tools/ubcg) follow the instructions for installing the program. 
    1. Create a folder named ReCo---- $mkdir ReCo
    2. save the downloaded UBCG program here
    3. UNZIP the folder: unzip UBCG_v3.zip
2. JAVA installation (Check for: jdk-8u5-linux-x64.tar.gz ):
    1. Unistalling JAVA completely: https://askubuntu.com/questions/84483/how-to-completely-uninstall-java#185250
    2. downloading: https://www.oracle.com/java/technologies/javase/javase8-archive-downloads.html (Check for: jdk-8u5-linux-x64.tar.gz )
    3. Re-Installing JAVA8 (JRE8): https://stackoverflow.com/questions/55920389/e-package-oracle-java8-installer-has-no-installation-candidate 
        note: use sudo in line of command
3. Install the dependencies using the sudo command. 
    1. sudo apt-get update -y
    2. sudo apt-get install -y prodigal
    3. sudo apt-get update -y
    4. sudo apt-get install -y mafft
    5. sudo apt-get update -y
    6. sudo apt-get install -y phipack
    7. sudo apt-get update -y
    8. sudo apt-get install -y raxml
    9. sudo apt-get update -y
    10. sudo apt-get install -y fasttree
4. After you have unzipped UBCG, you will notice a text file named 'programPath' in the UBCG directory. Open that file and edit it with following and save it thereafter:
    1.prodigal=/usr/bin/prodigal
    2.hmmsearch=/usr/bin/hmmsearch
    3.mafft=/usr/bin/mafft
    4.fasttree=/usr/bin/fasttree
    5.raxml=/usr/bin/raxmlHPC-PTHREADS
  
5. Now download the ReCo.py program and keep it in the unzipped UBCG folder present in the directory ReCo. 
6. RUN:- python3 ReCo.py (after running it will ask for the folders/ dstination of the folder having the scaffold/contigs of the sequences. ) 
7. Results are saved in the directory named CoRec.
