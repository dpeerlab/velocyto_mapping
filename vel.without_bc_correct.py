from seqc.sequence.encodings import DNA3Bit

from seqc.sequence import encodings
from seqc.sequence import barcodes as bcs
import pandas as pd
import subprocess
import os
import sys
import glob

import pysam
from seqc import read_array
from seqc.sequence.encodings import DNA3Bit


# Read array
prefix = sys.argv[1]
#prefix = 'DAC_B528-Kate+'

no_cores = 20

# Copy over genomic annotations and the barcodes
# aws s3 cp s3://dp-lab-data/genomes/mm38_long_polya/annotations.gtf ~/annot/
# aws s3 cp s3://dp-lab-data/barcodes/ten_x_v2/flat/737K-august-2016.txt ~/annot/
home_dir = ""
# Files
#data_dir = home_dir + f'{prefix}/{prefix}'
ra = read_array.ReadArray.load(prefix + '.h5')
in_bam = prefix + '_Aligned.out.sorted.bam'
out_bam = prefix + '_Aligned.tagged.bam'

# Read counts file to restrict the list of cells to use
tenx_bcs_bit = pd.read_csv(prefix + '_dense.csv', index_col=0).index


# Function to iterate over multiple alignments
def iter_multialignments(reader):
    """yields tuples of all alignments for each fastq record"""
    fq = [next(reader)]
    for record in reader:
        if record.qname == fq[0].query_name:
            fq.append(record)
        else:
            yield tuple(fq)
            fq = [record]
        yield tuple(fq)


# Open reader
print('Adding tags to bam file...')
reader = pysam.AlignmentFile(in_bam, "rb")
encoder = encodings.DNA3Bit()
with pysam.AlignmentFile(out_bam, "wb", header=reader.header) as outf:
    for read in reader.fetch():
        read_str = str(read)
        barcode = read_str.split(":")[1]
        umi = read_str.split(":")[2]

        print(barcode, umi)

        #barcode = DNA3Bit.decode(barcode).decode()

        if not DNA3Bit.encode(barcode) in tenx_bcs_bit:
            continue

        #umi = DNA3Bit.decode(umi).decode()

        #barcode_bit = encoder.encode(barcode)
        #barcode_corrected = bcs.find_correct_barcode(barcode_bit,tenx_bcs_bit)[0]
        #barcode_corrected_str = DNA3Bit.decode(barcode_corrected).decode()

        read.tags = read.tags + [('CB', barcode), ('UB', umi)]
        #read.tags = read.tags + [('CB', barcode_corrected_str), ('UB', umi)]
            # Write to bam
        outf.write(read)

reader.close()

del ra
import gc
gc.collect()

# Remove bam file to save space
#os.remove(in_bam)

# Sort tagged file
print('Sorting bam file...')
sorted_file = prefix + '_Aligned.tagged.sorted.bam'
pysam.sort("-o", sorted_file, '--threads', str(no_cores), out_bam)
# Remove unsorted file
#os.remove(out_bam)

# Index
args = f'samtools index {sorted_file}'
with subprocess.Popen(args, shell=True) as p:
    out, err = p.communicate()


# RNA velocity
print('Computing RNA velocity...')
out_dir = sys.argv[2] #f'{prefix}/velocity'
args = f'velocyto run -b /data/peer/chanj3/SCPC_transformation/ref/737K-august-2016.txt \
-o {out_dir} -@ {no_cores} -v {sorted_file} /data/peer/chanj3/SCPC_transformation/ref/annotations.gtf'

with subprocess.Popen(args, shell=True) as p:
    out, err = p.communicate()
