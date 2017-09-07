## provide a merged vcf files. The script generates a new merged
## vcf file that passes the missingness filter.

import optparse

def filter_missingness(infile,outfile,missvalue):
	vcf_infile = open(infile,'r')
	lines = vcf_infile.readlines()
	vcf_infile.close()
	
	vcf_outfile = open(outfile,'w')
	for i in lines:
		if i.startswith("#"):#comment lines
			vcf_outfile.write(i)
		else:
			col = i.strip().split()
			isolate_total = len(col) - 9 #total number of isolates
			isolate_missing = col.count("./.")#counts total no. of missing data
			missingness = float(isolate_missing)/isolate_total#calculates missingness value
			print ("missingness is ... " + str(missingness))
			if missingness <= float(missvalue):#passes the filter
				print ("passed the filter")
				vcf_outfile.write(i)
			else:#missingness value is greater than cutoff
				print ("didn't pass the filter")
	vcf_outfile.close()	

def main(infile,outfile,unique):
	print ("running main function")
	print "missingness cutoff is ... " + str(missvalue)
	filter_missingness(infile,outfile,missvalue)
 					
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
			missingness filter" )
	parser.add_option("-m", dest = "missvalue", type = "float",
		help = "the missingness cutoff: from 0 to 1")
	
    	# load the inputs
	(options , args) = parser.parse_args()
    
    	# process the inputs
	infile = options.infile
	outfile = options.outfile
	missvalue = options.missvalue
	main(infile,outfile,missvalue)
