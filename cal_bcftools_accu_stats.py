#calculates precision, recall, npv values for each read depth value 
#by using an input file which stores chrom coord, depth, TP/FP. 

import optparse


def read_infile(infile):
	'''
	open an input file, read all lines, store and return it as a list
	'''
	input_file = open(infile,'r')
	lines = input_file.readlines()
	input_file.close()
	#return a list contains all lines	
	return lines	

def get_max_depth(lines):
	'''
	return the maximum depth for all snps(TP/FP) from .accuracy.txt
	'''
	#use a list to store depth value
	tmp= []

	for i in lines:
		col = i.strip().split("\t")
		if i.startswith("#"):
			continue #skip the header
		else:
			dp = int(col[2])
			tmp.append(dp)
	max_depth = max(tmp)
	return max_depth


def filter_freq(lines,dp_cutoff):
	'''
	Extract all snps with TP/FP/TN/FN info
	with depth > dp_cutoff. Returns a list
	containing TP/FP/TN/FN info.
	'''
	#use a list to store TP/FP/TN/FN info
	list = []
	for i in lines:
		col = i.strip().split("\t")
		if i.startswith("#"):
			continue #skip the header
		elif (int(col[2]) < dp_cutoff+1):
			continue
		else:
			list.append(col[3])
	return list

def cal_pos_stats(lines,max_depth,pos_file,freq_file):
	'''
	Calculates precision and recall values
	'''
	outPosFile = open(pos_file,'w')
	outPosFile.write("MinReadCov\tPrecision\n") #write the header
	outFreqFile = open(freq_file,'w')
	outFreqFile.write("MinReadCov\tTP\tFP\n")
	for i in range(max_depth):
		info_list = filter_freq(lines,i)
		TP = info_list.count("TP")
		FP = info_list.count("FP")
		if (TP+FP != 0):
			outFreqFile.write(str(i+1) + "\t" + str(TP) + "\t" + str(FP) + "\n")
			precision = float(TP)/(TP+FP)#calculates precision value
			wl = str(i+1) + "\t" +  str(precision) + "\n"
			outPosFile.write(wl)
	outFreqFile.close()
	outPosFile.close()

def cal_neg_stats(lines,max_depth,neg_file):
        '''
        Calculates negative predictive value
        '''
        outfile = open(neg_file,'w')
        outfile.write("TotalReads\tNPV\n") #write the header
        for i in range(max_depth):
                info_list = filter_freq(lines,i)
                TN = info_list.count("TN")
                FN = info_list.count("FN")
                npv = float(TN)/(TN+FN)#calculates npv
                wl = str(i+1) + "\t" +  str(npv) + "\n"
                outfile.write(wl)
        outfile.close()

def main(infile,pos_file,neg_file,freq_file):
	lines = read_infile(infile)
	max_depth = get_max_depth(lines)
	cal_pos_stats(lines,max_depth,pos_file,freq_file)
	cal_neg_stats(lines,max_depth,neg_file)	
	
if __name__== '__main__':
	'''
	'''
    	# parser object for managing input options
	parser = optparse.OptionParser()

    	# essential data
	parser.add_option( '-i' , dest = 'infile' ,
                default = '' ,
                help = 'the input nanopolish_accuracy info file, which contains \
		TP, FP, FN' )
	parser.add_option( '-p' , dest = 'pos_file' ,
               	default = 'nanopolish_accu_stats_pos.txt' ,
	        help = 'specify output filename, which stores precision, recall values' )
	parser.add_option( '-n' , dest = 'neg_file' ,
                default = 'nanopolish_accu_stats_neg.txt' ,
                help = 'specify output filename, which stores negative predictive values' )
	parser.add_option( '-f' , dest = 'freq_file' ,
                default = '' ,
                help = 'specify output filename, which stores TP/FP frequency values' )

    	# load the inputs
	(options , args) = parser.parse_args()
    
    	# process the inputs
    	# note: all commandline inputs are str by default
	infile = options.infile
	pos_file = options.pos_file
	neg_file = options.neg_file
	freq_file = options.freq_file
	main(infile,pos_file,neg_file,freq_file)	
