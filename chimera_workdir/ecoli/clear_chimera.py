import sys

def write_fasta(fasta, seqs, names, width=80):
    with open(fasta, "w") as f:
        for idx in range(len(seqs)):
            seq = seqs[idx]
            name = names[idx]
            f.write(">{}\n".format(name))
            if width > 0:
                n = len(seq)
                for i in range(0, n, width):
                    f.write("{}\n".format(seq[i:min(i+width,n)]))
            else:
                f.write("{}\n".format(seq))

def read_fasta(fasta):
    seqs = []
    names = []
    with open(fasta, "r") as f:
        seqlines = []
        name = ""
        for line in f.readlines():
            if line.startswith(">"):
                if len(seqlines) > 0:
                    seqs.append("".join(seqlines))
                    names.append(name)
                    seqlines = []
                name = line.lstrip(">").rstrip()
            else:
                seqlines.append(line.rstrip())
        if len(seqlines) > 0:
            seqs.append("".join(seqlines))
            names.append(name)
    assert len(seqs) == len(names)
    return seqs, names


if __name__ == "__main__":

    if len(sys.argv) != 2:
        sys.stderr.write("usage: {} <reads.fa>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    fasta_fname = sys.argv[1]
    seqs, names = read_fasta(fasta_fname)

    for i in range(len(seqs)):
        if not "Rc" in names[i]:
            sys.stdout.write(">{}\n".format(names[i]))
            n = len(seqs[i])
            for j in range(0, n, 80):
                sys.stdout.write("{}\n".format(seqs[i][j:min(j+80,n)]))
    sys.stdout.flush()

