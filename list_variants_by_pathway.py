import vcf
import pprint
import subprocess
# more info: https://github.com/jamescasbon/PyVCF/

import argparse, textwrap, os

pp = pprint.PrettyPrinter(indent = 2)

parser = argparse.ArgumentParser(description='List variants by pathway')

parser.add_argument('input_vcf', metavar='i', nargs=1, help='input VCF (full path)')

args = parser.parse_args()

input_vcf = args.input_vcf[0]

# 1. build dictionary pathways-genes
d = {}
l = []
f = open('/path_to/h.all.v6.0.symbols.gmt')
lines = f.readlines()
for line in lines:
	v = line.split()
	pathway = v[0]
	genes = v[2:]
	d[pathway] = genes
	l.append(pathway)


f = open('/path_to/DRG.list')
genes = f.readlines()
d['DRG'] = [x.rstrip() for x in genes]
l.append('DRG')

f = open('/path_to/BROCA.list')
genes = f.readlines()
d['BROCA'] = [x.rstrip() for x in genes]
l.append('BROCA')

f.close()


# 2. iterate input VCF line by line and check if variant is in any pathway
dict_pathway_to_variant = {}
for pathway in l:
	dict_pathway_to_variant[pathway] = []

vcf_reader = vcf.Reader(open(input_vcf, 'r'))
for record in vcf_reader:
	current_gene = str(record.INFO['CSQ'][0]).split('|')[3]

	for pathway in l:
		if current_gene in d[pathway]:

			s = [str(record.CHROM), str(record.POS)]

			if record.ID == None: s.append('.')
			else: s.append(record.ID)

			s = s + [record.REF, str(record.ALT[0])]
			s = "\t".join(s)
			s = "^" + s

			sub1 = subprocess.Popen([ 'grep','-P', s, input_vcf ], stdout = subprocess.PIPE, shell = False)
			sub2 = subprocess.Popen( ['cut', '-f1-8'], stdin = sub1.stdout, stdout = subprocess.PIPE, shell = False)
			sub1.stdout.close()
			vcf_line = sub2.communicate()[0].rstrip()			
			num_carriers = record.num_het + record.num_hom_alt
			v = vcf_line.split('\t')
			vcf_line = v[0] + ':' + v[1] + ':' + v[3] + '>' + v[4] + '\t' + v[2] + '\t' + str(record.INFO['CSQ'][0]).split('|')[1] + '\t' + current_gene + '\t' + str(num_carriers)

			# vcf_line = subprocess.check_output(['grep','-P', s, input_vcf]).rstrip()

			# add current variant to dictionary
			dict_pathway_to_variant[pathway].append(vcf_line)


for pathway in l:
	print "/////////////////////"
	print pathway + ": " + str(len(dict_pathway_to_variant[pathway]))
	print "/////////////////////"
	print "Position\trsID\tVariant type\tGene\tNum carriers"
	for variant in dict_pathway_to_variant[pathway]:
		print variant


# python list_variants_by_pathway.py <input.vcf>


