"""
Remove multiallic SNPs, keep SNPs where REF isn't present in samps but SNP has two ALT alleles,
discard any biallelic SNPs with REF=N and one ALT.
This uses pandas chunks so that large tablefiles can be passed to script.

### usage
# python remove_multiallelic-keep_noREF.py infile outfile
###

### positional arguments
# infile = /path/to/file.txt output from GATK VariantsToTable
###

### assumes
# that output from VariantsToTable included CHROM, POS, AF, REF, ALT, and genotype fields from each sample
# that there is a non-zero number of non-multiallelic SNPs in VariantsToTable.txt infile
###
"""


import sys, pandas as pd, numpy as np, math
from coadaptree import uni
from os import path as op
from collections import Counter


def table(lst):
    """
    Count each item in a list.
    
    Returns:
    - Counter() with key = item, val = count
    """
    c = Counter()
    for x in lst:
        c[x] += 1
    return c


def get_chunks(tablefile):
    "Read in tab-delimited text file, return dataframe generator"
    basename = op.basename(tablefile)
    chunks = pd.read_csv(tablefile, sep='\t', chunksize=10000)
    return chunks, basename


def adjust_freqs(smalldf, alts):
    """
    For loci with REF not present, set ALT freqs with respect to the second ALT allele.
    
    Positional arguments:
    smalldf - pandas.dataframe; df where REF is not of the the two alleles present in samps
    
    Returns:
    df - smalldf with adjusted freqs in zeroth row
    """

    alt1, alt2 = alts
    df.loc[0, 'ALT'] = f"{alt1}+{alt2}"
    df.loc[0, 'AF'] = smalldf.loc[1, 'AF']
    return df 


def rm_multiallelic(tablefile):
    """
    Count CHROM-POS (locus) and keep only those with one ALT, discard if REF=N.
    
    Positional arguments:
    df - pandas.dataframe; currently filtered VariantsToTable output
    tf - basename of path to VariantsToTable output

    Returns:
    df - pandas.dataframe; non-multiallelic-filtered VariantsToTable output
    """
    print(f'removing multiallelic sites from {tablefile}')
    chunks, basename = get_chunks(tablefile)
    
    dfs = []
    for chunk in chunks:
        # give SNPs IDs by CHROM+POS
        chunk['locus'] = ["%s-%s" % (chrom,pos) for (chrom,pos) in zip(chunk["CHROM"], chunk["POS"])]
        # count SNPs, if count > 1 then it's a SNP with multiple ALT alleles
        loccount = table(chunk['locus'])
        # identify loci with exactly one ALT allele
        goodloci = [locus for locus in loccount if loccount[locus] == 1]
        # filter chunk for multiallelic (multiple lines), REF != N
        chunk = chunk[chunk['locus'].isin(goodloci)].copy()
        chunk = chunk[chunk['REF'] != 'N'].copy()
        # keep whatever is leftover after filtering
        dfs.append(chunk)
    # combine filtered chunks
    df = pd.concat(dfs)
    print(f'\t{tf} has {len(df.index)} good SNPs (non-multiallelic)')
    return df


def get_noref_snps(tablefile):
    """
    Isolate polymorphisms where REF isn't present in samps but SNP has two ALT alleles.
    
    Positional arguments:
    df - pandas.dataframe; current filtered VariantsToTable output
    
    Returns:
    dfs - list of loci (pandas.dataframes) to keep
    """
    print(f'identifying noREF SNPs from {tablefile} ...')
    chunks, basename = get_chunks(tablefile)

    dfs = []
    for chunk in chunks:
        # give SNPs IDs by CHROM+POS
        chunk['locus'] = ["%s-%s" % (chrom,pos) for (chrom,pos) in zip(chunk["CHROM"], chunk["POS"])]
        # count SNPs, if count > 1 then it's a SNP with multiple ALT alleles
        ncount = table(chunk['locus'])
        # identify loci with exactly two ALT alleles
        nloci = [locus for locus in ncount if ncount[locus] == 2]
        # reduce dataframe to loci with exactly two ALT alleles
        ndf = chunk[chunk['locus'].isin(nloci)].copy()
        # see which loci might have zero samps with REF allele
        for locus in nloci:
            smalldf = ndf[ndf['locus'] == locus].copy()
            if len(smalldf.index) == 2:  # redundant
                smalldf.index = range(len(smalldf.index))
                ref = smalldf.loc[0, 'REF']
                alts = smalldf['ALT'].tolist()
                keep = True
                gtcols = [col for col in smalldf.columns if '.GT' in col]
                for row in smalldf.index:
                    # get a string of alleles, delete "N" in case REF=N
                    gts = ''.join(smalldf.loc[row, gtcols].str.replace("/", "").tolist()).replace("N", "")
                    # if the REF is in the string, then it's a true multiallelic site and we discard
                    if ref in gts:
                        keep = False
                        break
                # if it seems to be a true multiallelic site and one of the ALTs is not *
                if keep is True and '*' not in alts:
                    newsmalldf = adjust_freqs(smalldf.copy(), alts)
                    dfs.append(pd.DataFrame(newsmalldf.loc[0,:]).T)
    print(f"\tfound {len(dfs)} SNPs where REF is not an allele in samps: {tf}")
    return dfs


def recombine(dfs, df):
    """
    Combine noREF SNPs with otherwise non-multallic SNPs.
    If zero noREF SNPs exist, return df as-is (which could have SNPs or not)
    """
    print('combining noREF SNPs and non-multiallelic SNPs ...')
    # if any noREF SNPs with two ALT alleles were found
    if len(dfs) > 0:
        # if tablefile (df) has normal SNPs (non noREF SNPs), combine all:
        if len(df.index) > 0:
            dfs.append(df)
            df = pd.concat(dfs)
        # else overwrite df with the noREF SNPs
        else:
            # if no normal SNPs, but > 0 noREF SNPs:
            df = pd.concat(dfs)
    # else: return df as-is
    return df

def main(tablefile, outfile):
    
    # isolate loci where REF is not of the the two alleles present in samps
    dfs = get_noref_snps(tablefile)
    
    # remove otherwise multiallelic sites
    df = rm_multiallelic(tablefile)
    
    # combine noREF SNPs with otherwise non-multiallelic SNPs
    df = recombine(dfs, df)
    
    # write file
    print(f'writing {len(df.index)} biallelic SNPs from {op.basename(tablefile)} to file: {outfile}')
    df.to_csv(outfile, sep='\t', index=False)
    

if __name__ == '__main__':
    thisfile, tablefile, outfile = sys.argv

    main(tablefile)
