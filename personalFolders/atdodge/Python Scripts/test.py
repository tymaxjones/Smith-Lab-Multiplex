import primer3
from Bio import SeqIO
import re

file_in = '../Sequence Files/test_input.fasta'
regex_match = r'[A]{5,}|[T]{5,}|[C]{5,}|[T]{5,}'
for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
	print(seq_record.id)
	primer_count = 0
	for i in range(0, len(seq_record.seq) - 17):
		current_slice = str(seq_record.seq[i:i + 18])
		if re.search(regex_match, current_slice):
			pass
		else:
			tm = primer3.calcTm(current_slice)
			position = "{}-{}".format(i+1, i+18)
			if tm > 55 and tm < 65 and current_slice[-1] != "T":
				print(current_slice + "\t{}\t{}".format(position,"Tm = " + str(tm)))
				primer_count += 1
print("Sequence length:", len(seq_record.seq))
print("Possible primers:", primer_count)
print()
