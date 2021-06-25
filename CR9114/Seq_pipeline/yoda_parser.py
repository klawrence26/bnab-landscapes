
import sys
from collections import OrderedDict


import regex

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def rc(seq):
	return "".join(complements.get(base, base) for base in reversed(seq))

def OR(xs):
	return "(" + "|".join(["(?:"+x+")" for x in xs]) + ")"

# to correct multiplexing indices
def make_corrector(options, num_mismatch=0):
	checkers = [regex.compile("("+o+"){s<=" + str(num_mismatch) + "}") for o in options]
	def corrector(match):
		current_error_count = 10000
		current_best_i = ""
		for (i,c) in enumerate(checkers):
			m = c.fullmatch(match)
			if m:
				error_count = m.fuzzy_counts[0] + m.fuzzy_counts[1] + m.fuzzy_counts[2]
				if(error_count == 0):
					return i+1
				elif (error_count < current_error_count and error_count <= num_mismatch):
					current_best_i = i+1
					current_error_count = error_count
		if (current_best_i != ""):
			return current_best_i
		else:
			return 0
	return corrector

P5_index = ['TAGATCGC',
	'CTCTCTAT',
	'TATCCTCT',
	'AGAGTAGA',
	'GTAAGGAG',
	'ACTGCATA',
	'AAGGAGTA',
	'CTAAGCCT',
	'GGCTACTC',
	'CCTCAGAC',
	'TCCTTACG',
	'ACGCGTGG',
	'GGAACTCC',
	'TGGCCATG',
	'GAGAGATT',
	'CGCGGTTA',
	'GACCGCCA',
	'TAAGATGG',
	'ATTGACAT',
	'AGCCAACT',
	'TACTAGGT',
	'TCACGGTT',
	'TGTAATGA',
	'CACGTCAG',
	'CTGAATTC',
	'CGTACCGG',
	'GATGACGG',
	'TATAGACG',
	'GTCATTGA',
	'GCATCGTT',
	'AGGTTGAC',
	'TGAAACTG',
	'CAATCACA',
	'ACATGCAA',
	'ATCGCGCC',
	'TCGGTTAA']

correct_P5 = make_corrector(P5_index)

P5indexDict = {idx:(i+1) for i,idx in enumerate(P5_index)}

yoda_col_index = [
	"TTGATCG",
	'XXXXXXXXXXXXXXXXXXXXXXXXX', 
	'XXXXXXXXXXXXXXXXXXXXXXXXX',
	'XXXXXXXXXXXXXXXXXXXXXXXXX',
	"GGTATCA",
	'XXXXXXXXXXXXXXXXXXXXXXXXX',
	'XXXXXXXXXXXXXXXXXXXXXXXXX',
	"GACGGAACTC",
	"ATCGCAGGCAT",
	"TCACGACTAGTA"]


correct_yoda_col = make_corrector(yoda_col_index,2)

yoda_row_index = [
	"TACCTGA",
	"AGAACCAT",
	'XXXXXXXXXXXXXXXXXXXXXXXXX',
	'XXXXXXXXXXXXXXXXXXXXXXXXX',
	"TCGGTGGTACG",
	"GTTCAACCGATT",
	"CGTACAA",
	'XXXXXXXXXXXXXXXXXXXXXXXXX']

correct_yoda_row = make_corrector(yoda_row_index,2)

P7_index = ['TCGCCTTA',
	'CTAGTACG',
	'TTCTGCCT',
	'GCTCAGGA',
	'AGGAGTCC',
	'CATGCCTA',
	'GTAGAGAG',
	'CCTCTCTG',
	'AGCGTAGC',
	'CAGCCTCG',
	'TGCCTCTT',
	'TCCTCTAC',
	'ATTACAAT',
	'GAATGATC',
	'CGATCGGT',
	'AATAACGG',
	'TAGAAGAA',
	'GTCAGGTA',
	'GCGGTCCT',
	'AATCGGAC',
	'AACTCGTG',
	'GGCCGTGG',
	'TTACATGT',
	'AGTTAACA']
	

correct_P7 = make_corrector(P7_index)

P7indexDict = {idx:(i+1) for i,idx in enumerate(P7_index)}



# Constant regions
molecular_barcode = "(.{8})"

forward_9114 = '(?e)(?:AAAGCCTCCGGTGGTACC){e<=2}'
const1 = '(?e)(?:TACGCCATAAGCTGGGTAAGGCAGGCCCCGGGGCAGGGCTTGGAGTGGATGGGCGGCATC){s<=3}'
const2 = '(?e)(?:CCCATATTCGGC){s<=1}'
const3 = '(?e)(?:TATGCACAGAAATTTCA){s<=2}'

const4rc = '(?e)(?:GTCCGC){s<=1}' 
const5rc = '(?e)(?:CAGTTCCATGTAAGCTGT){s<=3}'
const6rc = '(?e)(?:CAGACT){s<=1}'
const7rc = '(?e)(?:GTAGACTGCTGTATCTTCAGA){s<=2}'
const8rc = '(?e)(?:ATAGTAGTAGTAGTTACCGTGTCGTGCACA){s<=2}'
reverse_9114rc = '(?e)(?:GACCCCAGACATCCATACC){e<=2}'

mut = '(.{3})'



yoda_R1 = ( molecular_barcode + 
	    OR(yoda_col_index) + "{s<=2}" +
		 forward_9114 +
		 mut + mut + mut +
		 const1 +
		 mut +
		 const2 +
		 mut + mut + mut +
		 const3 )
		 

yoda_R2 = ( molecular_barcode +
		OR(yoda_row_index) + "{s<=2}" +
		reverse_9114rc  +
		mut +
		const8rc +
		mut +
		const7rc +
		mut +
		const6rc +
		mut +
		const5rc +
		mut + mut + mut + mut +
		const4rc +
		mut )
		 

yoda_R1re = regex.compile(yoda_R1)
yoda_R2re = regex.compile(yoda_R2)



nreads = 0
n_match = 0
r1_discarded = 0
r2_discarded = 0
for line in sys.stdin:
	nreads += 1
	output=OrderedDict([
		    (x,"") for x in ["illumina_index_P5","illumina_index_P7","col","row","genotype","mb1","mb2"]
	    ])
	
	splitline = line.split("	")
	if len(splitline) < 4:
		continue

	R1 = splitline[2]
	R2 = splitline[3]
 
	m1 = yoda_R1re.match(R1)
	if m1:
		m2 = yoda_R2re.match(R2)
		
		
		if m2:
			
			n_match += 1
			output["illumina_index_P7"] = P7indexDict[rc(splitline[0])]
			output["illumina_index_P5"] = P5indexDict[splitline[1]]
			output["mb1"] = m1.groups()[0]
			output["col"] = correct_yoda_col(m1.groups()[1])
			
			
			output["mb2"] = m2.groups()[0]
			output["row"] = correct_yoda_row(m2.groups()[1])
			geno_fwd = m1.groups()[2]+'.'+m1.groups()[3]+'.'+m1.groups()[4]+'.'+m1.groups()[5]+'.'+m1.groups()[6]+'.'+m1.groups()[7]+'.'+m1.groups()[8]
			geno_rev = rc(m2.groups()[2]+'.'+m2.groups()[3]+'.'+m2.groups()[4]+'.'+m2.groups()[5]+'.'+m2.groups()[6]+'.'+m2.groups()[7]+'.'+m2.groups()[8]+'.'+m2.groups()[9]+'.'+m2.groups()[10])
			output["genotype"] = geno_fwd+'.'+geno_rev
			
					
			print(*output.values(), sep="	", file=sys.stdout)	

		else:
			r2_discarded += 1
			#print(splitline[3])
			continue
	else:
		r1_discarded += 1
		#print(splitline[2])
		continue
 
	
 
print("Discarded " + str(r1_discarded) + " + " + str(r2_discarded) + " out of " + str(nreads) + " reads.", file=sys.stderr)
print("Parsed "+str(n_match/nreads*100)+"% of reads",file=sys.stderr)
print(str(n_match)+'\t'+str(nreads)+'\t'+str(n_match/nreads*100),file=sys.stdout)

