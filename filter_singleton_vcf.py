## This script removes singletons from a merged vcf file
## it only accepts homozygous SNPs (0/0, 1/1)
import optparse
import re

def filter_singleton(infile,outfile):
	vcf_infile = open(infile,'r')
	lines = vcf_infile.readlines()
	vcf_infile.close()
	
	vcf_outfile = open(outfile,'w')
	for i in lines:
		if i.startswith("#"):#comment lines
			vcf_outfile.write(i)
		else:
			singleton_counts = re.findall("1/1",i)
			if len(singleton_counts) != 1:#if it's not singleton
				vcf_outfile.write(i)
	vcf_outfile.close()	

def main(infile,outfile):
	print("running main function")
	filter_singleton(infile,outfile)
 					
if __name__== '__main__':
	'''
	main function
	'''
    	# parser object for managing input options
	parser = optparse.OptionParser()

    	# essential data
	parser.add_option( "-i", dest = "infile",
		default = "",
		help = 'provide input .vcf file' )
	
	parser.add_option( "-o", dest = "outfile",
                default = "",
                help = "name the output .vcf file which passes the \
			singleton filter" )
	
    	# load the inputs
	(options , args) = parser.parse_args()
    
    	# process the inputs
    	# note: all commandline inputs are str by default
	infile = options.infile
	outfile = options.outfile
	main(infile,outfile)
