import vcf
# more info: https://github.com/jamescasbon/PyVCF/

import argparse, textwrap, os


# parser = argparse.ArgumentParser(description='Filter variants for: PTMs 95% (LOFTEE LoF=HC), missense CADD PHRED > 20, ExAC MAF NFE < 1%, repeat regions and monomorphic loci.')
parser = argparse.ArgumentParser(
	formatter_class = argparse.RawDescriptionHelpFormatter,
	description = textwrap.dedent('''\
		Retains the following variants:
			- Tier1: PTMs 95% (LOFTEE LoF=HC)
			- Tier2: missense CADD PHRED > 20
			Both Tier1 and Tier2 must pass the following QC:
			- ExAC MAF NFE < 1%
			- removes variants falling in repeat regions (FILTER=repeat_region)
			- removes monomorphic loci.
			'''))

parser.add_argument('input_vcf', metavar='i', nargs=1, help='input VCF (full path)')
parser.add_argument('--g', metavar = 'GENE', help = 'output variants only in this gene')

# parser.add_argument('output_vcf', metavar='o', nargs=1, help='output VCF (full path)')

# parser.add_argument('integers', metavar='N', type=int, nargs='+',
# 					help='an integer for the accumulator')

# parser.add_argument('--sum', dest='accumulate', action='store_const',
# 					const=sum, default=max,
# 					help='sum the integers (default: find the max)')

# print args.accumulate(args.integers)

args = parser.parse_args()

input_vcf = args.input_vcf[0]
if os.path.dirname(input_vcf) == '':
	full_path = "."
else:
	full_path = os.path.dirname(input_vcf)


g = args.g

if g is None:
	output_vcf_Tier1 = full_path + "/" + os.path.basename(input_vcf)[:-3] + "Tier1.vcf"
	output_vcf_Tier1_Tier2 = full_path + "/" + os.path.basename(input_vcf)[:-3] + "Tier1_Tier2.vcf"
else:
	output_vcf_Tier1 = full_path + "/" + os.path.basename(input_vcf)[:-3] + g + ".Tier1.vcf"
	output_vcf_Tier1_Tier2 = full_path + "/" + os.path.basename(input_vcf)[:-3] + g + ".Tier1_Tier2.vcf"


vcf_reader = vcf.Reader(open(input_vcf, 'r'))
vcf_writer_Tier1 = vcf.Writer(open(output_vcf_Tier1, 'w'), vcf_reader)
vcf_writer_Tier1_Tier2 = vcf.Writer(open(output_vcf_Tier1_Tier2, 'w'), vcf_reader)

for record in vcf_reader:

	if str(record.ALT[0]) == '*':
		continue
	
	mut_type = str(record.INFO['CSQ'][0]).split('|')[1] # missense
	impact = str(record.INFO['CSQ'][0]).split('|')[2] # HIGH, MODERATE
	gene = str(record.INFO['CSQ'][0]).split('|')[3]
	if "CLNSIG" not in record.INFO:
		clin_sig = "-"
	else:
		clin_sig = str(record.INFO['CLNSIG'][0]) # Likely_benign

	LoF = str(record.INFO['CSQ'][0]).split('|')[27] # LoF=HC/LC, only for PTMs
	position = str(record.INFO['CSQ'][0]).split('|')[30] # LoF_info=POSITION:0.0751633986928105
	ExAC_AF_NFE = str(record.INFO['CSQ'][0]).split('|')[40]
	CADD_PHRED = str(record.INFO['CSQ'][0]).split('|')[44]
	REVEL = str(record.INFO['CSQ'][0]).split('|')[46]
	gnomAD_exomes_AF_nfe = str(record.INFO['CSQ'][0]).split('|')[48]

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
	l_hom_alts = ",".join([h.sample for h in record.get_hom_alts()])


	is_Tier1 = False
	is_Tier2 = False

	if impact == 'HIGH' and LoF == "HC":
		is_Tier1 = True
	elif impact == 'MODERATE' and mut_type.find('missense') != -1 and (CADD_PHRED == "" or float(CADD_PHRED) > 20):
		is_Tier2 = True

	# check also if variant is monomorphic
	n = len(record.samples)
	is_monomorphic = (record.num_het == n or record.num_hom_ref == n or record.num_hom_alt == n)
	FILTER = ",".join(record.FILTER) # if falls in a repeat region, then FILTER='repeat_region'
	pass_filter = (is_Tier1 or is_Tier2) and (ExAC_AF_NFE == "" or float(ExAC_AF_NFE) <= 0.01) and FILTER == "" and not is_monomorphic
	if not g is None:
		pass_filter = pass_filter and (gene == g)

	if pass_filter:
		# print record.CHROM + '\t' + str(record.POS) + '\t' + str(record.ID) + '\t' + record.REF + '\t' + str(record.ALT[0]) + '\tFILTER=' + FILTER + '\tTYPE=' + mut_type + '\tIMPACT=' + impact + '\tGENE=' + gene + '\tCLIN_SIG=' + clin_sig + '\tLoF=' + LoF + '\tPosition=' + position + '\tExAC_MAF=' + ExAC_AF_NFE + '\tgnomAD_exomes_AF_nfe=' + gnomAD_exomes_AF_nfe + '\tCADD_PHRED=' + CADD_PHRED + '\tREVEL=' + REVEL + "\tHET=" + l_het + "\tHOM_ALT=" + l_hom_alts
		if is_Tier1:
			vcf_writer_Tier1.write_record(record)
			vcf_writer_Tier1_Tier2.write_record(record)
		elif is_Tier2:
			vcf_writer_Tier1_Tier2.write_record(record)


# python filter_variants.py input.vcf --g ATM




