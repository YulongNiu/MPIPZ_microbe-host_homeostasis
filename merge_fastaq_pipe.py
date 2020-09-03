# simple merger for 2 fastq files
# written for python 2.7
#
# please make sure input files are sorted
# and have the same length
# and same file format

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import gzip
import sys

parser = argparse.ArgumentParser("merges 2 fastq files into one by concatenating the sequences and qualitysequences\n make sure both files have the same format!")
parser.add_argument("fileA", type=str, help="fastq file A")
parser.add_argument("fileB", type=str, help="fastq file B")
parser.add_argument("--out", default=None, help="name outputfile - otherwise output to std-out. If file name ends with .gz file will be compressed automatically.")
#parser.add_argument("--gzip-input-files", default=False, action="store_true", help="set flag if files are compressed")
parser.add_argument("--format", default="fastq", type=str, help="specify your fastq format - default is fastq")


args = parser.parse_args()
fq_format = args.format

# read files

try:
	if args.fileA.endswith(".gz"): #automatic gz scanning   args.gzip_input_files:
		handleA = gzip.open(args.fileA,"rt")
		handleB = gzip.open(args.fileB,"rt")

		fileA = SeqIO.parse(handleA,fq_format)
		fileB = SeqIO.parse(handleB,fq_format)
	else:
		fileA = SeqIO.parse(args.fileA,fq_format)
		fileB = SeqIO.parse(args.fileB,fq_format)

except Exception as e:
	print("something went wrong reading files \n"+e.message)
	sys.exit(0)

# merge and output

if args.out is None:
	output = sys.stdout
elif args.out.endswith(".gz"):
	output = gzip.open(args.out,"wt")
else:
	output = open(args.out,"w")

while True:
	try:
		recA = next(fileA)
		recB = next(fileB)
		if recA.id != recB.id:
			sys.stderr.write("found non matching ids, skipping, make sure files are sorted\n")
			continue

		out_record = SeqRecord(recA.seq+recB.seq, id = recA.id, name = recA.name, description = recA.description+"|"+recB.description)
		out_record.letter_annotations["phred_quality"] = recA.letter_annotations["phred_quality"]+recB.letter_annotations["phred_quality"]
	
		output.write(out_record.format(fq_format))

	except StopIteration as e:
		break

	except Exception as e:
		sys.stderr.write("Something went wrong.\n"+e.message)
		sys.exit()


