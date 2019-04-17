import vcf
# more info: https://github.com/jamescasbon/PyVCF/

import argparse, textwrap, os

# parser = argparse.ArgumentParser(description='Filter variants for: PTMs 95% (LOFTEE LoF=HC), missense CADD PHRED > 20, ExAC MAF NFE < 1%, repeat regions and monomorphic loci.')
parser = argparse.ArgumentParser(
	formatter_class = argparse.RawDescriptionHelpFormatter,
	description = textwrap.dedent('''\
		Write a txt file with the variants filtered by the optional arguments:
			- Tier1: PTMs 95% (LOFTEE LoF=HC)
			- Tier2: missense CADD PHRED > 20
		All variants are filtered by the following QC:
			- ExAC MAF NFE < 1%
			- removes variants falling in repeat regions (FILTER=repeat_region)
			- removes monomorphic loci
			'''))

parser.add_argument('input_vcf', help = 'input VCF (full path)', metavar = '<input_vcf>', nargs = 1)
parser.add_argument('-g', help = 'output variants only in this gene', metavar = 'GENE')
parser.add_argument('-t1', help = 'include Tier1 variants', action = 'store_true')
parser.add_argument('-t2', help = 'include Tier2 variants', action = 'store_true')
parser.add_argument('-vcf', help = 'Additionally, output in VCF format', action = 'store_true')

# parser.add_argument('output_vcf', metavar='o', nargs=1, help='output VCF (full path)')

# parser.add_argument('integers', metavar='N', type=int, nargs='+',
# 					help='an integer for the accumulator')

# parser.add_argument('--sum', dest='accumulate', action='store_const',
# 					const=sum, default=max,
# 					help='sum the integers (default: find the max)')

# print args.accumulate(args.integers)

args = parser.parse_args()

input_vcf = args.input_vcf[0]
t1 = args.t1
t2 = args.t2
g = args.g
output_in_vcf_format = args.vcf


if os.path.dirname(input_vcf) == '':
	full_path = "."
else:
	full_path = os.path.dirname(input_vcf)

if os.path.basename(input_vcf).endswith('vcf.gz'):
	output_vcf = full_path + "/" + os.path.basename(input_vcf)[:-6] + "filtered.vcf"
	output_txt = full_path + "/" + os.path.basename(input_vcf)[:-6] + "filtered.txt"
else:
	output_vcf = full_path + "/" + os.path.basename(input_vcf)[:-3] + "filtered.vcf"
	output_txt = full_path + "/" + os.path.basename(input_vcf)[:-3] + "filtered.txt"


vcf_reader = vcf.Reader(open(input_vcf, 'r'))

output_txt = open(output_txt, 'w')
output_txt.write("CHROM\tPOS\trsID\tREF\tALT\tFILTER\tTYPE\tIMPACT\tGENE\tCLIN_SIG\tLoF\tLoF_info\tExAC_MAF_NFE\tgnomAD_exomes_NFE\tCADD_PHRED\tREVEL\tHET_CARRIERS\tHOM_ALT_CARRIERS\n")


if output_in_vcf_format:
	vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)

for record in vcf_reader:


	if str(record.ALT[0]) == '*':
		continue
	
	mut_type = str(record.INFO['CSQ'][0]).split('|')[1]
	impact = str(record.INFO['CSQ'][0]).split('|')[2]
	gene = str(record.INFO['CSQ'][0]).split('|')[3]
	clin_sig = "-" if "CLNSIG" not in record.INFO else str(record.INFO['CLNSIG'][0])
	LoF = str(record.INFO['CSQ'][0]).split('|')[27] if str(record.INFO['CSQ'][0]).split('|')[27] != "" else "-"
	LoF_info = str(record.INFO['CSQ'][0]).split('|')[30] if str(record.INFO['CSQ'][0]).split('|')[30] != "" else "-"
	ExAC_AF_NFE = str(record.INFO['CSQ'][0]).split('|')[40] if str(record.INFO['CSQ'][0]).split('|')[40] != "" else "-"
	CADD_PHRED = str(record.INFO['CSQ'][0]).split('|')[44] if str(record.INFO['CSQ'][0]).split('|')[44] != "" else "-"
	REVEL = str(record.INFO['CSQ'][0]).split('|')[46] if str(record.INFO['CSQ'][0]).split('|')[46] != "" else "-"
	gnomAD_exomes_AF_nfe = str(record.INFO['CSQ'][0]).split('|')[48] if str(record.INFO['CSQ'][0]).split('|')[48] != "" else "-"

	# print record.get_hets()
	# print record.get_hom_alts()
	# print record.samples
	# print record.num_het
	# print record.num_hom_alt
	# print record.num_hom_ref
	# print record.is_snp
	# print record.is_indel
	# print record.is_monomorphic

	l_het = ",".join([h.sample for h in record.get_hets()])
	if l_het == "": l_het = "-"
	l_hom_alts = ",".join([h.sample for h in record.get_hom_alts()])
	if l_hom_alts == "": l_hom_alts = "-"

	is_Tier1 = False
	is_Tier2 = False

	if impact == 'HIGH' and LoF == "HC":
		is_Tier1 = True
	elif impact == 'MODERATE' and mut_type.find('missense') != -1 and (CADD_PHRED == "-" or float(CADD_PHRED) > 20):
		is_Tier2 = True

	# check monomorphic
	n = len(record.samples)
	is_monomorphic = (record.num_het == n or record.num_hom_ref == n or record.num_hom_alt == n)
	# check repeated region
	FILTER = ",".join(record.FILTER)

	pass_filter = (ExAC_AF_NFE == "-" or float(ExAC_AF_NFE) <= 0.01) and FILTER == "" and not is_monomorphic

	# if str(record.POS) == "9471381":
	# 	print record.INFO
	# 	print is_Tier1, LoF_info
	# 	print t1
	# 	exit()

	if not g is None:
		pass_filter = pass_filter and (gene == g)

	if pass_filter:
		if (t1 and is_Tier1) or (t2 and is_Tier2) or (not t1 and not t2):
			# print record.CHROM + '\t' + str(record.POS) + '\t' + str(record.ID) + '\t' + record.REF + '\t' + str(record.ALT[0]) + '\tFILTER=' + FILTER + '\tTYPE=' + mut_type + '\tIMPACT=' + impact + '\tGENE=' + gene + '\tCLIN_SIG=' + clin_sig + '\tLoF=' + LoF + '\tLoF_info=' + LoF_info + '\tExAC_MAF=' + ExAC_AF_NFE + '\tgnomAD_exomes_AF_nfe=' + gnomAD_exomes_AF_nfe + '\tCADD_PHRED=' + CADD_PHRED + '\tREVEL=' + REVEL + "\tHET=" + l_het + "\tHOM_ALT=" + l_hom_alts
			output_txt.write(record.CHROM + '\t' + str(record.POS) + '\t' + str(record.ID) + '\t' + record.REF + '\t' + str(record.ALT[0]) + '\t' + FILTER + '\t' + mut_type + '\t' + impact + '\t' + gene + '\t' + clin_sig + '\t' + LoF + '\t' + LoF_info + '\t' + ExAC_AF_NFE + '\t' + gnomAD_exomes_AF_nfe + '\t' + CADD_PHRED + '\t' + REVEL + '\t' + l_het + '\t' + l_hom_alts + '\n')
			if output_in_vcf_format:
				vcf_writer.write_record(record)

if output_in_vcf_format:
	vcf_writer.close()
output_txt.close()


# python filter_variants.py -h
# python filter_variants.py chr2.vcf.gz -t1



