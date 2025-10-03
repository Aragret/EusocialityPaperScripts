# --------------------------------------------------------------------------------------------
# ------------- Main script for the CAFE5 analysis of the Eusociality project ----------------
# --------------------------------------------------------------------------------------------

# ------------- Required data ----------------------------------------------------------------

# The following input files are required: 
# tcal_blattodea_tree.tre
# gene_hog_small_count_eusoc.txt
# gene_hog_large_count_eusoc.txt

# ------------- Making of folder architecture ------------------------------------------------

mkdir CAFE
cd CAFE
mkdir Runs_for_convergence/
mkdir Runs_for_convergence/gamma/
mkdir Runs_for_convergence/cafe_optimal/
mkdir trees_for_optimum/
mkdir slurmK
mkdir slurmM
mkdir lists

# ------------- Search optimal gamma parameter -----------------------------------------------

sbatch --array=1-500 Job-gammasearch_lambdafixed.sh

# ------------- Extracting result of optimal gamma parameter --------------------------------

echo -e "rate\tlnL\tlb1\talp" > gamma_selection.txt
ls slurmK > slurmK_list.txt
for fi in $(cat slurmK_list.txt); do
	rate=$(echo $(grep "Command line:" slurmK/$fi)| cut -d "/" -f8)
	lnL=$(echo $(grep "Score (-lnL):" slurmK/$fi | tail -1) | cut -d " " -f3)
	lb1=$(echo $(grep "Best match" slurmK/$fi | tail -1)| cut -d " " -f4 | cut -d "," -f1)
	alp=$(echo $(grep "Best matches are:" slurmK/$fi)| cut -d " " -f4 | cut -d "," -f2)
	echo -e "${rate}\t${lnL}\t${lb1}\t${alp}" >> gamma_selection.txt
done

# ------------- Result visualisation: gamma --------------------------------------------------

# The results of the gamma selection were visualised using the R script:
# 1.r_script_for_gamma_selection.R

# ------------- Preparing multiple tree ------------------------------------------------------

cp tcal_blattodea_tree.tre tree_M0.tre

# The blattodea_tree_M0.tre was made manually using tree_M0.tre

# for the models:
for i in {1..56}
do
	sed "s/1/2/${i}" tree_M0.tre > trees_for_optimum/tree_M${i}.tre
done
ls trees_for_optimum > "M_list.txt"

# The M_rep_list.txt file was made manually with excel and text editor 

# ------------- Search optimal lambda parameter ----------------------------------------------

sbatch --array=1-2800 Job-lambdasearch_gammafixed_repi.sh  # (50 repetition, 56x50)
     #=> slurmM/slurm files

# ------------- Extracting result of optimal lambda parameter --------------------------------

ls slurmM > slurmM_list.txt
echo -e "mod\tlnL\tlb1\tlb2" > Mlambda_selection.txt
for fi in $(cat slurmM_list.txt); do
	mod=$(echo $(grep "tree_M" slurmM/$fi)| cut -d "_" -f10 | cut -d "." -f1)
	lnL=$(echo $(grep "Score (-lnL):" slurmM/$fi | tail -1) | cut -d " " -f3)
	lb1=$(echo $(grep "Best matches are:" slurmM/$fi)| cut -d " " -f4 | cut -d "," -f1)
	lb2=$(echo $(grep "Best matches are:" slurmM/$fi)| cut -d " " -f4 | cut -d "," -f2)
	echo -e "${mod}\t${lnL}\t${lb1}\t${lb2}" >> Mlambda_selection.txt
done

# ------------- lambda clustering ------------------------------------------------------------

# The lambda per branch were clustered using the R script:
# 2.script for lambda selection by kmeans.R

# ------------- Preparing trees for model selection ------------------------------------------

cp tree_M0.tre tree_Mod.tre
for i in {1..56}; do
	sed -i "s/:1/:M${i}/1" tree_Mod.tre
done
for u in {2..6}; do
	cp tree_Mod.tre optimal_tree_l${u}.tre

	j=1
	for k in $(cat clusters_l${u}.txt); do
		sed -i "s/:M${j}/:${k}/1" optimal_tree_l${u}.tre
		let "j+=1"
	done
done


# The l_k_rep.txt file was made manually using Excel and text editor
# example: 
#1/1/1
#2/1/1
#...
#6/3/1
#2/1/2
#...
#6/3/50


# ------------- Model selection ---------------------------------------------------------------

sbatch --array=1-900 Job-opt_cafe_l_k_rep.sh                                     #RUNNING 50 rep

# ------------ Extracting result of model selection -------------------------------------------


ls slurmO > slurmO_list.txt
echo -e "mod\tno_rate\tKi\tno_cluster\tLi\tlbi\tlnL\tlambda" > optimal_selection.txt

for fi in $(cat slurmO_list.txt); do
	mod=$(echo $(grep "Command line:" slurmO/$fi| cut -d "/" -f8))
	i=$(echo $mod | cut -d "l" -f2)
	Ki=$(echo $mod | cut -d "_" -f1)
	no_rate=$(echo $mod | cut -d "_" -f1 | cut -d "k" -f2)
	lnL=$(echo $(grep "Score (-lnL):" slurmO/$fi | tail -1) | cut -d " " -f3)
	canc=$(echo $(grep "CANCELLED" slurmO/$fi))
	if [ -z "${canc}" ]; then for ((u=1;u<=i;u++)); do lb=$(echo $(grep "Best match" slurmO/$fi | tail -1)| cut -d " " -f4 | cut -d "," -f${u}); echo -e "${mod}\t${no_rate}\t${Ki}\t${i}\tL${i}\tl${u}\t${lnL}\t${lb}" >> optimal_selection.txt; done; fi
done

# ------------- Result visualisation: optimal selection ----------------------------------------

# The results of the optimal selection were visualised using the R script:
# 3.make_the_plot_for_selection.R

# ------------ Selection -----------------------------------------------------------------------

# We selected K3L1. 
# Note: the model with K=2 and L=2 was tested further but did not converge
# We are now repeating the model for K3 L1 50 times with the set value of lambdas.
# Then use the script to rune the large gene families

sbatch --array=1-50 Job-k3_L1_best.sh

#Note:
#When using gamma models, there is an indication on what are the gene families that have more 
#than 20% failure during optimisation. That means when CAFE is looking for a value of lambda
#and alpha, it may run into a dead end and tries another value. It will compare the values it 
#found and get the lowest set of values based on likelihood. But the best match are always a 
#set of values that did not fail. Wether they are to take depend on the convergence of these 
#values over multiple tests. If the values converge, then they are satisfying.



