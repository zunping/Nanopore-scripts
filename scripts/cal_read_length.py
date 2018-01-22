import optparse

'''
calculate read length from a given fasta file
'''

def cal_read_length(fasta,out_filename):
	infile = open(fasta,'r')
	lines = infile.readlines()
	infile.close()
	
	outfile = open(out_filename,'w')
	for i in lines:
		if not i.startswith(">"):
			length = len(i.strip())
			outfile.write(str(length)+"\n")
	outfile.close()

if __name__== '__main__':
    	# parser object for managing input options
	parser = optparse.OptionParser()

    	# essential data
	parser.add_option( '-f' , dest = 'fasta' ,
		default = '' ,
		help = 'the input fasta file' )

	parser.add_option( '-o' , dest = 'out_filename' ,
                default = '' ,
                help = 'the output filename' )

    	# load the inputs
	(options , args) = parser.parse_args()
    
    	# process the inputs
    	# note: all commandline inputs are str by default
	fasta = options.fasta
	out_filename = options.out_filename
	cal_read_length(fasta,out_filename)
