from re import template
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

my_seq = Seq("GACTACTTGGAATGA")

for index, letter in enumerate(my_seq):
  print("%i %s" % (index, letter))
print(my_seq[0]) #first letter
print(my_seq[-1]) #last letter

print(my_seq.count("A"))

# Purin asoslari DNK da necha foizligini bilish uchun:
print("================================")
foiz = 100 * float(my_seq.count("A")+my_seq.count("T")/len(my_seq))
print(foiz)
print(GC123(my_seq))

## Switch string
print("================================")
seq_str = str(my_seq)
print(seq_str)

fasta_format_str = ">Name\n%s\n" % seq_str
print(fasta_format_str)

## Slicing
print("================================")
print(my_seq[0:13])
print(my_seq[0::3])

## Concatenating or adding sequences
print("================================")
protein_seq = Seq("EVRNAK")
dna_seq = Seq("ACGT")
print(protein_seq + dna_seq)

list_of_seqs = [Seq("ACGT"), Seq("AACC"), Seq("GGTT")]
concatenated = Seq("")
for s in list_of_seqs:
  concatenated += s
print(concatenated)

contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
spacer = Seq("N"*10)
print(spacer.join(contigs))

## Changing case
print("================================")
dna_seq_two = Seq("agctAGCA")
print("CTA" in dna_seq_two)
print("CTA" in dna_seq_two.upper())

## Nucleotide sequences and (reverse) complements
print("================================")
my_seq_three = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
print(my_seq_three.complement())
print(my_seq_three.reverse_complement())
print(my_seq_three.complement_rna())
print(my_seq_three[::-1])

## Transcription
print("================================")
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print(coding_dna)

templated_dna = coding_dna.reverse_complement()
print(templated_dna)

messenger_rna = coding_dna.transcribe()
print(coding_dna.complement())
print(messenger_rna)

messenger_rna_two = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
# 5'-3' ga qaytaradi
print(messenger_rna_two.back_transcribe())

## Translation
print("================================")
# oqsil = messenger_rna_two.translate()
oqsil = coding_dna.translate(table="Vertebrate Mitochondrial")
print(oqsil)

oqsil_two = coding_dna.translate(table=2)
print(oqsil_two)

print(coding_dna.translate(to_stop=True))
print(coding_dna.translate(table=2, to_stop=True))

gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA"
          "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT"
          "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT"
          "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT"
          "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")

print(gene.translate(table="Bacterial"))
print(gene.translate(table="Bacterial", to_stop=True))
print(gene.translate(table="Bacterial", cds=True))

##  Translation Tables
print("================================")
from Bio.Data import CodonTable
# standart_table = CodonTable.unambiguous_dna_by_name["Standart"]
# mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

standard_table = CodonTable.unambiguous_dna_by_id[1]
mito_table = CodonTable.unambiguous_dna_by_id[2]
other_table = CodonTable.unambiguous_dna_by_id[3]
print(standard_table)
print(mito_table)
print(other_table)

print(mito_table.stop_codons)
print(mito_table.start_codons)
print(mito_table.forward_table["ACG"])

## Comparing Seq objects
print("================================")
seq = Seq("ACT")
print("ACT" == seq)
print(seq == "ACT")

## Sequences with unknown sequence contents
print("================================")
uncnown_seq = Seq(None, 10)
print(len(uncnown_seq))
# print(uncnown_seq)

## MutableSeq objects
"""
Oddiy Python qatori singari, Seq ob'ekti "faqat o'qiladi" yoki Python terminologiyasida o'zgarmasdir.
Seq ob'ektining satr kabi ishlashini xohlashdan tashqari, bu ham foydali standartdir, chunki ko'plab 
biologik ilovalarda siz ketma-ketlik ma'lumotlarini o'zgartirmasligingizga ishonch hosil qilishni xohlaysiz:
"""
# my_seq[5] = "T"

"""
However, you can convert it into a mutable sequence (a MutableSeq object) and do pretty much anything
you want with it:
"""
from Bio.Seq import MutableSeq
mutable_seq = MutableSeq(my_seq)
print(mutable_seq)
test = mutable_seq[2]="T"
print(mutable_seq)

"""
Shu bilan bir qatorda, MutableSeq ob'ektini to'g'ridan-to'g'ri satrdan yaratishingiz mumkin:
"""
mutable_seq_two = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
mutable_seq_two[5] = "C"
print(mutable_seq_two)