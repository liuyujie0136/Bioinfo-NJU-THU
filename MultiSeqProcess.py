f = open("MultiSeq_p53.txt", 'r')
seq_dict = {}
seq_name = None
for line in f:
    if line[0] == '>':
        seq_name = line[1:].strip()
        seq_dict[seq_name] = ''
    else:
        if seq_name:
            seq_dict[seq_name] = seq_dict[seq_name] + line.strip()
print(seq_dict)
f.close()
