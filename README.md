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
    a. Create a folder named ReCo---- $mkdir ReCo
    b. save the downloaded UBCG program here
    c. UNZIP the folder: unzip UBCG_v3.zip
2. JAVA installation (Check for: jdk-8u5-linux-x64.tar.gz ):
  -Unistalling JAVA completely: https://askubuntu.com/questions/84483/how-to-completely-uninstall-java#185250
  -Re-Installing JAVA8 (JRE8): https://stackoverflow.com/questions/55920389/e-package-oracle-java8-installer-has-no-installation-candidate (Check for: jdk-8u5-linux-x64.tar.gz )
3. Install the dependencies using the sudo command. 
  sudo apt-get update -y
  sudo apt-get install -y prodigal
  sudo apt-get update -y
  sudo apt-get install -y mafft
  sudo apt-get update -y
  sudo apt-get install -y phipack
  sudo apt-get update -y
  sudo apt-get install -y raxml
  sudo apt-get update -y
  sudo apt-get install -y fasttree
4. After you have unzipped UBCG, you will notice a text file named 'programPath' in the UBCG directory. Open that file and edit it with following and save it thereafter:
  prodigal=/usr/bin/prodigal
  hmmsearch=/usr/bin/hmmsearch
  mafft=/usr/bin/mafft
  fasttree=/usr/bin/fasttree
  raxml=/usr/bin/raxmlHPC-PTHREADS
  
5. Now download the ReCo.py program and keep it in the unzipped UBCG folder present in the directory ReCo. RUN:- python3 ReCo.py (after running it will ask for the folders/ dstination of the folder having the scaffold/contigs of the sequences. ) 
6. Results are saved in the directory named CoRec.
