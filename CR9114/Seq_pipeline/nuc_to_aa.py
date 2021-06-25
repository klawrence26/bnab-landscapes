
import sys
import regex


aa_dict_list = [{'TTT':'F','AGT':'S'},
					{'TCT':'S','AAC':'N'},
					{'TCC':'S','AAC':'N'},
					{'ATC':'I','AGT':'S'},
					{'ACC':'T','TCT':'S'},
					{'GCG':'A','ACA':'T'},
					{'AAT':'N','GCT':'A'},
					
					{'ACC':'T','TCT':'S'},
					{'AAG':'K','ATA':'I'},
					{'AGC':'S','TTT':'F'},
					{'ACT':'T','TCC':'S'},
					{'TCC':'S','ATT':'N'},
					{'AGC':'S','AAT':'N'},
					{'AGA':'R','ACT':'T'},
					{'TAC':'Y','TTT':'F'},
					{'TAC':'Y','AGC':'S'}]

binary_dict_list = [{'TTT':'0','AGT':'1'},
					{'TCT':'0','AAC':'1'},
					{'TCC':'0','AAC':'1'},
					{'ATC':'0','AGT':'1'},
					{'ACC':'0','TCT':'1'},
					{'GCG':'0','ACA':'1'},
					{'AAT':'0','GCT':'1'},
					
					{'ACC':'0','TCT':'1'},
					{'AAG':'0','ATA':'1'},
					{'AGC':'0','TTT':'1'},
					{'ACT':'0','TCC':'1'},
					{'TCC':'0','AAT':'1'},
					{'AGC':'0','AAT':'1'},
					{'AGA':'0','ACT':'1'},
					{'TAC':'0','TTT':'1'},
					{'TAC':'0','AGC':'1'}]
					
germline_codons = ['TTT','TCT','TCC','ATC','ACC','GCG','AAT',				
						'ACC','AAG','AGC','ACT','TCC','AGC','AGA','TAC','TAC']
somatic_codons = ['AGT','AAC','AAC','AGT','TCT','ACA','GCT',
					'TCT','ATA','TTT','TCC','AAT','AAT','ACT','TTT','AGC']						

num_muts = 16

# loop through parsed reads to find errors
n_without_errors = 0
n_with_uncorrected_errors = 0
n_with_corrected_errors = 0
for line in sys.stdin:
	corrected_flag = False
	error_flag = False
	splitline = line.split("	")
	p5,p7,col,row,geno,mb1,mb2 = splitline
	
	mut_list = geno.split('.')
	
	geno_to_print = []
	
	for i in range(num_muts):
		if mut_list[i] in binary_dict_list[i].keys():
			geno_to_print.append(binary_dict_list[i][mut_list[i]])
				
		else:
			geno_to_print.append('X')
			error_flag = True

	if error_flag: n_with_uncorrected_errors += 1
	elif corrected_flag: n_with_corrected_errors += 1
	else: 
		n_without_errors += 1	
		print(p5,p7,col,row,''.join(geno_to_print),mb1,mb2,sep='\t', file=sys.stdout,end='')
	
#print('Fraction of reads with errors: ',float(n_with_uncorrected_errors)/(float(n_with_corrected_errors)+float(n_with_uncorrected_errors)+float(n_without_errors)),file=sys.stderr)

#print(str(n_without_errors)+'\t'+str(n_with_corrected_errors)+'\t'+str(n_with_uncorrected_errors)+'\t'+str(float(n_with_uncorrected_errors)/(float(n_with_corrected_errors)+float(n_with_uncorrected_errors)+float(n_without_errors))),file=sys.stdout)
