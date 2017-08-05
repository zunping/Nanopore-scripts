import optparse

def filter_missingness(infile,outfile,missvalue):
	vcf_infile = open(infile,'r')
	lines = vcf_infile.readlines()
	vcf_infile.close()
	
	vcf_outfile = open(outfile,'w')
	for i in lines:
		if i.startswith("#"):
			vcf_outfile.write(i)
		else:
			col = i.strip().split()
		########idx = col.index("GT:AD:DP:GQ:PL")
			isolate_total = len(col) - 9 #total number of isolates
			isolate_missing = col.count("./.")
			missingness = float(isolate_missing)/isolate_total
			print "missingness is ... " + str(missingness) 
			if missingness <= float(missvalue):
				print "pass filter"
				vcf_outfile.write(i)
			else:
				print "fail filter"
	vcf_outfile.close()	

def main(infile,outfile,unique):
	print "running main function"
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
                help = "name the output .vcf file which pass the \
			missingnes filter" )
	parser.add_option("-m", dest = "missvalue", 
		help = "the missingness cutoff: from 0 to 1")
	
    	# load the inputs
	(options , args) = parser.parse_args()
    
    	# process the inputs
    	# note: all commandline inputs are str by default
	infile = options.infile
	outfile = options.outfile
	missvalue = options.missvalue
	main(infile,outfile,missvalue)
