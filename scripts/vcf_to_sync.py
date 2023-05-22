import argparse
from cyvcf2 import VCF

from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""description""")
parser.add_argument('--input', required=True, dest='Input', type=str, help="Input file in vcf format")
parser.add_argument('--output', required=True, dest='Output', type=str, help="Output file in sync format.")

args = parser.parse_args()


#InputFile=open(args.Input,"r")
OutputFile=open(args.Output,"w")

base_dict = {'A':0, 'T':1, 'C':2, 'G':3}

###
# define conversion function

def convert_to_sync(ref_allele, var_allele, ref_count, alt_count):
    ct_ref=ref_count
    ct_var=alt_count
    sync_line = "0:0:0:0:0:0".split(":")
    sync_line[base_dict.get(ref_allele)]=ct_ref
    sync_line[base_dict.get(var_allele)]=ct_var
    return ":".join(sync_line)

##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 20">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">

for variant in VCF(args.Input): # or VCF('some.bcf')
    ref_var = variant.REF
    alt_var = variant.ALT # e.g. REF='A', ALT=['C', 'T']

    # make first 3 columns:
    col_chrom = variant.CHROM
    col_pos = variant.end

    # get depths for all samples:
    rd = variant.format('RD')
    ad = variant.format('AD')
    sync_samps = [None] * len(rd)
    for samp in range(len(rd)):
    	sync_samps[samp] = convert_to_sync(ref_var, alt_var[0], ' '.join(map(str, rd[samp])), ' '.join(map(str, ad[samp])))

    # make new complete line
    line_out = (col_chrom + "\t" + str(col_pos) + "\t" +  ref_var + "\t" + '\t'.join(map(str,sync_samps)))

    OutputFile.write(line_out+'\n')


OutputFile.close()
