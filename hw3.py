from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse
from collections import Counter


class Kmer_spectres:

    def __init__(self):
        self.file = ''
        self.sequence = ''
        self.kmer_dict = {}
        self.counter = 0
        self.spectres = None
        self.value = []
        self.freq = []


    def kmer_search(self, file, kmer_size, quality):
        for record in SeqIO.parse(file, 'fastq'):
            num = 0
            qual = record.letter_annotations['phred_quality']
            for i in qual:
                if i < quality:
                   num = 1
            if num == 1:
                continue
            local_seq = record.seq
            for index in range(len(local_seq) - kmer_size + 1):
                kmer = local_seq[index:(index + kmer_size)]
                if kmer in self.kmer_dict:
                    self.kmer_dict[kmer] = self.kmer_dict[kmer] + 1
                else:
                    self.kmer_dict[kmer] = 1

    def spectres_plot(self):
        self.spectres = Counter(self.kmer_dict.values()).most_common()
        for i in self.spectres:
            self.value.append(i[0])
            self.freq.append(i[1])
        plt.bar(self.value, self.freq, width=5, color='red')
        #print(self.spectres)
        self.max = 0
        for i in range(len(self.spectres)):
            if self.spectres[i][1] > self.max:
                self.max = self.spectres[i][1]
        plt.ylim((0, self.max))
        plt.tight_layout()
        plt.show()


    def spectres_thresh(self, thresh):
        fig, ax = plt.subplots()
        ax.bar(self.value, self.freq, width=5, color='red')
        plt.ylim((0, self.max))
        plt.axvline(x=int(thresh))
        plt.tight_layout()
        plt.show()

    def genome_size(self, thresh):
        count = 0
        count1 = 0
        count2 = 0
        for i in range(len(self.spectres)):
            if self.spectres[i][0] > int(thresh):
                count += self.spectres[i][0] * self.spectres[i][1]
                count1 += self.spectres[i][0]
                count2 += 1
        count1 = count1 / count2
        size = count / count1
        print("___________Genome size:___________\n", int(size))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Kmer spectres')
    parser.add_argument('-f', '--file', help='input file in the fastq format', type=str, required=True)
    parser.add_argument('-s', '--size', help='kmer size', type=int, default=15)
    parser.add_argument('-q', '--quality', help='quality threshold', type=int, default=20)
    parser.add_argument('-n', '--noise', help = 'Threshold of noise value', type = int, default=5)

    args = parser.parse_args()
    file = args.file
    size = args.size
    qual = args.quality
    noise_treshold = args.noise

    data = Kmer_spectres()
    data.kmer_search(file, size, qual)
    data.spectres_plot()
    data.spectres_thresh(noise_treshold)
    data.genome_size(noise_treshold)