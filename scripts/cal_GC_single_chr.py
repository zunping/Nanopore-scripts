import optparse
from Bio import SeqIO
import re


'''
calculate the GC content for a given 
chr/scaffold from a given fasta file
'''


def create_dict(fasta):
	dict = {}
	for i in SeqIO.parse(fasta,"fasta"):
		seq_name = i.id
		sequence = i.seq
		dict[seq_name]=sequence
	return dict

def main(fasta,chr,window,outfile):
	with open(outfile,'w') as f:
		print "Now creaing dictionary for chr..." + chr
		dict = create_dict(fasta)
		print "dictionary is ready"
		sequence = dict[chr]
		seq_len = len(sequence)
		index = seq_len/window
		for i in range(index):
			window_seq = sequence[i*window:(i+1)*window]
			gc = cal_GC(window_seq)	
			writeline = str(i*window) + "\t" + gc + "\n"
			f.write(writeline)
		last_seq = sequence[index*window:]
		gc = cal_GC(last_seq)
		writeline = str(index*window) + "\t" + gc + "\n"
		f.write(writeline)

'''
returns GC content if the seq isn't full of "N", otherwise
returns zero.
'''
def cal_GC(seq):
	gcCount = len(re.findall("[GCgc]", str(seq)))
	totalBaseCount = len(re.findall("[GCTAgcta]", str(seq)))
	if totalBaseCount != 0:
		GC_content = 100*(float(gcCount)/totalBaseCount)
	else:
		GC_content = 0
	return str(GC_content) 
	
if __name__== '__main__':
	'''
	calculate the GC content for a chr/scaffold from a given fasta file
	'''
    	# parser object for managing input options
	parser = optparse.OptionParser()

    	# essential data
	parser.add_option( '-f', dest = 'fasta',
		default = '',
		help = 'the input fasta file' )
	
	parser.add_option( '-c' , dest = 'chr',
                default = '' ,
                help = 'the chromosome name' )

	parser.add_option( '-w' , dest = 'window', type="int",
                default = 100 ,
                help = 'the sliding window size. default is 100bp' )

	parser.add_option( '-o' , dest = 'outfile' ,
                default = '' ,
                help = 'the output file containing GC content for each sequence in a give fasta file' )
	
    	# load the inputs
	(options , args) = parser.parse_args()
    
    	# process the inputs
    	# note: all commandline inputs are str by default
	fasta = options.fasta
	chr = options.chr
	window = options.window
	outfile = options.outfile
	main(fasta,chr,window,outfile)
