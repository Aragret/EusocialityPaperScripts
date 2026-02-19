#!/usr/bin/python3

########## multi fasta file library##################

##contains functions: read_fasta, write_fasta, rev_comp, translate, hydrophobicity, Tm_est, GC_content

################################################################################################################### 

### read multi-fasta file and return dictionary containing name:sequence pairs
def read_fasta(file_name):
	with open(file_name, "r") as f:  #read file into f
		seq_dict = {} #create empty sequence dictionary, we want names as keys and sequences as values
		seq = "" #create empty string for appending sequence lines. Important if sequence is broken over several lines.
		for line in f:
			if line[0] == ">":
				if bool(seq): #if seq has been assigned and we are now at a "name" line,
					     #first assign previous sequence to previous name
					seq_dict[name] = seq
					seq = "" # then empty the seq again for the next sequence
				name = line.strip()
				name = name.split()[0][1:] #get name from first word in line which begins with ">"
			else:
				seq = seq + line.strip() #keep adding sequence until another name line is read
		seq_dict[name] = seq
	return seq_dict


################################################################################################################### 

### write multi-sequence fasta using dictionary of name:sequence pairs
def write_fasta(seq_dict, file_name):
	with open(file_name, "w") as f:
		for i in seq_dict:
			f.write(">" + str(i) + "\n" + str(seq_dict[i]) + "\n")

################################################################################################################### 

### calculate reverse complement of multiple sequences from dictionary
def rev_comp(seq_dict):
	rc_dict = {}
	for i in seq_dict:
		seq_dict[i] = seq_dict[i].upper()
		rc_list = []
		for base in seq_dict[i]:
			if base == "A":
				rc_list.insert(0,"T")
			elif base == "T":
				rc_list.insert(0,"A")
			elif base == "C":
				rc_list.insert(0,"G")
			elif base == "G":
				rc_list.insert(0,"C")
		rc_dict[i] = "".join(rc_list)
	return rc_dict

################################################################################################################### 

### Translate DNA sequence in dictionary to AA seqs
def translate(seq_dict):
	codon_dict = {	"TTT":"F", "TCT":"S", "TAT":"Y", "TGT":"C",
			"TTC":"F", "TCC":"S", "TAC":"Y", "TGC":"C",
			"TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*",
			"TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W",
			"CTT":"L", "CCT":"P", "CAT":"H", "CGT":"R",
			"CTC":"L", "CCC":"P", "CAC":"H", "CGC":"R",
			"CTA":"L", "CCA":"P", "CAA":"Q", "CGA":"R",
			"CTG":"L", "CCG":"P", "CAG":"Q", "CGG":"R",
			"ATT":"I", "ACT":"T", "AAT":"N", "AGT":"S",
			"ATC":"I", "ACC":"T", "AAC":"N", "AGC":"S",
			"ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R",
			"ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R",
			"GTT":"V", "GCT":"A", "GAT":"D", "GGT":"G",
			"GTC":"V", "GCC":"A", "GAC":"D", "GGC":"G",
			"GTA":"V", "GCA":"A", "GAA":"E", "GGA":"G",
			"GTG":"V", "GCG":"A", "GAG":"E", "GGG":"G"  }
	prot_dict = {}
	for seq in seq_dict:
		prot_seq = []
		seq_dict[seq] = seq_dict[seq].upper()
#		if seq_dict[seq][0:3] != "ATG":
#			print(seq , "is not a coding sequence, must begin with start codon. ")
#			break
		if len(seq_dict[seq])%3 != 0:
			print(seq , "is not complete, sequence must be divisibile by 3.")
			break
		else:
			for i in range(0,len(seq_dict[seq])-3,3):
				prot_seq.append(codon_dict[seq_dict[seq][i:i+3]])
			prot_dict[seq] = "".join(prot_seq)
	return prot_dict

###################################################################################################################

### Calculate hydrophobicity

def hydrophobicity(prot_seq,window = 21): # a protein sequence; optional size of sliding window
	hydro_dict = {	"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,
			"M":1.9,"A":1.8,"G":-0.4,"T":-0.7,"S":-0.8,
			"W":-0.9,"Y":-1.3,"P":-1.6,"H":-3.2,"E":-3.5,
			"Q":-3.5,"D":-3.5,"N":-3.5,"K":-3.9,"R":-4.5}
	prot_seq = prot_seq.upper()
	hydro_list = []
	for i in range(len(prot_seq)-window):
		total = 0
		if prot_seq[i] not in hydro_dict:
			print("Sequence contains invalid amino acid.")
			break
		else:
			for aa in prot_seq[i:i+window]:
				total += float(hydro_dict[aa])
			hydro_list.append(total/window)

	return hydro_list


###################################################################################################################

### Tm estimator

def Tm_est(seq): #enter dna sequence
	seq = seq.upper()
	temp = None
	for base in seq:
		if base not in "ACGT":
			print("contains unallowed character:",base)
			break
	else:
		if len(seq) < 14:
			temp = 4 * (seq.count("G") + seq.count("C")) + 2*(seq.count("A") + seq.count("T"))
		else:
			temp = 64.9 + 41 * ((seq.count("G") + seq.count("C") - 16.4)/len(seq))
	return temp

################################################################################################################### 

### GC content calculator

def gc_calc(seq): #enter dna sequence
	seq = seq.upper()
	gc_count = None
	for base in seq:
		if base not in "ACGT":
			print("contains unallowed character:",base)
			break
	else:
		gc_count = (seq.count("G") + seq.count("C")) / len(seq)
	return gc_count
	
###################################################################################################################


if __name__ == "__main__":
	chrebp_dict = read_fasta("CHREBP.fasta")
#	rc_chrebp_dict = rev_comp(chrebp_dict)
#	write_fasta(rc_chrebp_dict,"CHREBP_rc.fasta")
#	chrebp_pep_dict = translate(chrebp_dict)
#	write_fasta(chrebp_pep_dict, "CHREBP_pep.fasta" )
#	chrebp_hydro = hydrophobicity(chrebp_pep_dict["Mnat_15238"])
#	with open("chrebp_hydro.txt","w") as f:
#		for i in chrebp_hydro:
#			f.write(str(i)+"\n")
	print("Tm is ", Tm_est(chrebp_dict["Mnat_15238"]))
	print("GC content is", gc_calc(chrebp_dict["Mnat_15238"]))
