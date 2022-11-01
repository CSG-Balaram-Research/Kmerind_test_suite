import pyfastx
import sys, os

#First command line param K
#Second relativepath to the file
def get_inputs():
    n_params = len(sys.argv)
    assert(n_params == 3)
    k = int(sys.argv[1]) 
    file = sys.argv[2]
    return k, file


#reference https://pypi.org/project/pyfastx/, https://claresloggett.github.io/python_workshops/improved_kmers.html
def get_counts(k, file):
    counts = {}
    fq = pyfastx.Fastx(file)
    for rec in fq:
        #rec = (name,seq,qual,comment)
        gnome = rec[1] 
        num_kmers = len(gnome) - k + 1
        for i in range(num_kmers):
            kmer = gnome[i:i+k]
            if kmer not in counts.keys():
                counts[kmer] = 1
            else:
                counts[kmer] += 1
    return counts


def write_out(counts, file):
    out_file = file + ".counts"
    if(os.path.exists(out_file)):
        os.remove(out_file)
    f = open(out_file, "w")
    for kmer, count in counts.items():
        ans = kmer + "@" + str(count) + "\n"
        f.write(ans)
    f.close()


# function simmulation c main function
def main():
    K, input_file = get_inputs()
    counts = get_counts(K, input_file)
    write_out(counts, input_file)

main()
