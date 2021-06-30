import primer3
from Bio import SeqIO
from tkinter import Tk
from tkinter.filedialog import askopenfilename


Tk().withdraw()
test_file = askopenfilename()
file_parsed = SeqIO.parse(open(test_file, mode='r'), 'fasta')

for seq_record in file_parsed:
	print(primer3.bindings.designPrimers(
		{
			'SEQUENCE_ID': seq_record.id,
			'SEQUENCE_TEMPLATE': str(seq_record.seq)
		},
		{
			'PRIMER_OPT_SIZE': 18,
			'PRIMER_MAX_SIZE': 25,
			'PRIMER_OPT_TM': 60.0
		}
	))

print(primer3.bindings.designPrimers(
	{
		'SEQUENCE_ID': 'TEST',
		'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG'
	},
	{
		'PRIMER_OPT_SIZE': 18,
		'PRIMER_MAX_SIZE': 25,
		'PRIMER_OPT_TM': 60.0,
		'PRIMER_NUM_RETURN': 7
	}
))
