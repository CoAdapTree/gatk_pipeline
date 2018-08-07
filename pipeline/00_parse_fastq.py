from __future__ import division
from Bio import SeqIO
import sys
import os
from os import path as op
import gzip
os.system('source $HOME/.bashrc')

thisfile,f,pooldir = sys.argv

def batch_iterator(iterator, batch_size):
    entry = True  # Make sure loops once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

            
record_iter = SeqIO.parse(gzip.open(f),"fastq")
batches = batch_iterator(record_iter, 10000)
for i, batch in enumerate(batches):
    print 'started'
    f = op.basename(f).replace("fastq.gz","")
    filename = op.join(pooldir,"%s_%s.fastq" % (f,
                                                 str(i).zfill(4)))
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fastq")
    print("Wrote %i records to %s" % (count, filename))

