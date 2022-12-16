# from Bio.Seq import Seq

# my_seq = Seq("AAGCTTGCAGCC")

# print(f"List: {my_seq}\n"
#       f"Complement RNA: {my_seq.complement_rna()}\n"
#       f"Complement: {my_seq.complement()}\n"
#       f"Reverse Complement: {my_seq.reverse_complement()}")

from Bio import SeqIO

## Simple FASTA parsing example
# for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
#   print(seq_record.id)
#   print(repr(seq_record.seq))
#   print(len(seq_record))

## Simple GenBank parsing example
for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
  print(seq_record.id)
  print(repr(seq_record.seq))
  print(len(seq_record))