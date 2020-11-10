from Bio import SeqIO
f = open("MultiSeq_p53.txt", "r")
for record in SeqIO.parse(f, "fasta"):
    print(record.id)
    print(record.name)
    print(record.seq)
f.close()
