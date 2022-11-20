import sys
import numpy as np
import random

def random_chromosome(length):
    return "".join(random.choice("ACGT") for i in range(length))

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def reverse_complement(s):
    return s.translate(comp_tab)[::-1]

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

def circular_slice(s, i, l):
    n = len(s)
    assert l < n
    i %= n
    startpos = i
    if i+l >= n:
        endpos = i+l-1
        cslice = s[i:endpos+1]
    else:
        endpos = i+l-n-1
        cslice = s[i:] + s[:endpos+1]
    return (cslice, startpos, endpos)

def create_subread(genome, startpos, readlen):
    endpos = startpos+readlen-1
    readseq = genome[startpos:endpos+1]
    startpos += 1
    endpos += 1
    if bool(random.randint(0,1)):
        readseq = reverse_complement(readseq)
        startpos, endpos = endpos, startpos
    return (readseq, startpos, endpos)


def create_read(genome, readlen, chimera):
    records = []
    splitpos = -1
    startpos = random.randint(0, len(genome)-readlen-1)
    if not chimera:
        records.append(create_subread(genome, startpos, readlen))
    else:
        records = []
        sub1len = random.randint(1, readlen-2)
        splitpos = sub1len+1
        records.append(create_subread(genome, startpos, sub1len))
        sub2len = readlen - sub1len
        startpos = random.randint(0, len(genome)-sub2len-1)
        records.append(create_subread(genome, startpos, sub2len))
    return records, splitpos

def main(ref_filename, reads_filename, depth, read_length_mean, read_length_sd, chimera_fraction):
    chroms, chrom_names = read_fasta(ref_filename)
    genome = chroms[0]
    genome_name = chrom_names[0].split()[0]
    genome_length = len(genome)
    num_reads = int((genome_length * depth) / read_length_mean)

    reads = []
    read_names = []

    for i in range(num_reads):
        readlen = 0
        while readlen <= 0:
            readlen = int(np.random.normal(read_length_mean, read_length_sd))
        records, splitpos = create_read(genome, readlen, bool(random.uniform(0,1) < chimera_fraction))
        if len(records) == 1:
            readseq, startpos, endpos = records[0]
            readname = "R{}\tlength={}\tcoords={}[{}..{}]".format(i+1, len(readseq), genome_name, startpos, endpos)
            reads.append(readseq)
            read_names.append(readname)
        else:
            readseq1, startpos1, endpos1 = records[0]
            readseq2, startpos2, endpos2 = records[1]
            readseq = readseq1 + readseq2
            readname = "Rc{},{}\tlength={}\tcoords1={}[{}..{}]\tcoords2={}[{}..{}]".format(i+1, splitpos, len(readseq), genome_name, startpos1, endpos1, genome_name, startpos2, endpos2)
            reads.append(readseq)
            read_names.append(readname)


    write_fasta(reads_filename, reads, read_names)

    return 0

def main_multiple_chroms(ref_filename, reads_filename, depth, read_length_mean, read_length_sd, chimera_fraction):
    chroms, chrom_names = read_fasta(ref_filename)
    num_chroms = len(chroms)
    genome_length = sum(len(chrom) for chrom in chroms)
    num_reads = int((genome_length * depth) / read_length_mean)

    reads = []
    read_names = []

    for i in range(num_reads):
        chrom_id = random.randint(0, num_chroms-1)
        chrom = chroms[chrom_id]
        chrom_name = chrom_names[chrom_id]
        chrom_len = len(chrom)
        readpos = random.randint(0, chrom_len-1)
        while True:
            readlen = int(np.random.normal(read_length_mean, read_length_sd))
            if readlen > 0:
                break
        readseq, startpos, endpos = circular_slice(chrom, readpos, readlen)
        readrev = bool(random.randint(0,1))
        if startpos < endpos:
            coords = "[{}..{}]".format(startpos, endpos)
        else:
            coords = "[{}..)++[..{}]".format(startpos, endpos)
        readname = "{}\tlength={}\tcoords={}\trev={}\tchrom={}".format(i+1, readlen, coords, readrev, chrom_name)
        if readrev:
            readseq = reverse_complement(readseq)
        reads.append(readseq)
        read_names.append(readname)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        exe = sys.argv[0]
        sys.stderr.write("usage: {} <ref.fa> <reads.fa> <depth> <read length mean> <read length deviation> <chimera fraction>\n".format(exe))
        sys.stderr.flush()
        sys.exit(-1)
    ref_filename = sys.argv[1]
    reads_filename = sys.argv[2]
    depth = int(sys.argv[3])
    read_length_mean = float(sys.argv[4])
    read_length_sd = float(sys.argv[5])
    chimera_fraction = float(sys.argv[6])

    retval = main(ref_filename, reads_filename, depth, read_length_mean, read_length_sd, chimera_fraction)
    sys.exit(retval)
