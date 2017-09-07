import optparse
#reads a merged vcf file, which contains SNPs info for multiple samples.
#for each sample, the script concatenates all SNPs and then 
#generates a fasta file that contains concatenanted sequence for each sample.


def create_dict(infile):
	in_file = open(infile,'r')
	lines = in_file.readlines()
	in_file.close()
	
	#creat a dict to store all variants
	#each key is a isolate name
	#each value is the concatenated SNPs for each isolate
	dict = {}
	for i in lines:
		if i.startswith("#CHROM"):
			isolate_names = i.strip().split("\t")#extract isolates' name
	#initiate dictionary
	for c in isolate_names[9:]:
		dict[c] = ""
	return dict,isolate_names


def concat_var(infile):
	dict = create_dict(infile)[0]
	isolate_names = create_dict(infile)[1]
	in_file = open(infile,'r')
	lines = in_file.readlines()
	in_file.close()
	
	for i in lines:
		if i.startswith("#"):
			continue

		col = i.strip().split("\t") #line contains snp information
		ref = col[3] #reference allele
		alt = col[4] #alternative allele
		length = len(col)	
		
		if len(ref) > 1 or len(alt) > 1:
			continue
		else:
			for c in range(9,length):
				isolate = isolate_names[c]
				snp = col[c] #extract snp info
				value = dict[isolate]
				if snp.startswith("0/0"):
					value += ref
				elif snp.startswith("1/1"):
					value += alt
				elif snp.startswith("./."):
					value += "-"
				dict[isolate] = value
	return dict


def print_dict(infile,outfile):
	dict = concat_var(infile)
	out_file = open(outfile,'w')
	for i in dict:
		final_key = ">" + i + "\n"
		seq = dict[i] + "\n"
		out_file.write(final_key)
		out_file.write(seq)
	out_file.close()

def main(infile,outfile):
	print "running main function"
	print_dict(infile,outfile)
 					
if __name__== '__main__':
	'''
	main function
	'''
    	# parser object for managing input options
	parser = optparse.OptionParser()

    	# essential data
	parser.add_option( "-i", dest = "infile",
		default = "",
		help = 'provide input .vcf file, heterozygous sites masked (./.) \
			no singleton' )
	
	parser.add_option( "-o", dest = "outfile",
                default = "",
                help = "name the output .fasta file which contains all variants" )
	
    	# load the inputs
	(options , args) = parser.parse_args()
    
    	# process the inputs
    	# note: all commandline inputs are str by default
	infile = options.infile
	outfile = options.outfile
	main(infile,outfile)
