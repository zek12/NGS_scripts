import vcf
# more info: https://github.com/jamescasbon/PyVCF/

import argparse, textwrap, os

# parser = argparse.ArgumentParser(description='Filter variants for: PTMs 95% (=LOFTEE HC), missense CADD PHRED > 30, ExAC MAF NFE < 1%')
parser = argparse.ArgumentParser(
	formatter_class = argparse.RawDescriptionHelpFormatter,
	description = textwrap.dedent('''\
		Filter variants for:
			- PTMs 95% (LOFTEE LoF=HC)
			- missense CADD PHRED > 30
			Also removes variants with ExAC MAF NFE < 1%, monomorphic sites, or FILTER column (repeat region)
			'''))

parser.add_argument('input_vcf', metavar='i', nargs=1, help='input VCF (full path)')

# parser.add_argument('output_vcf', metavar='o', nargs=1, help='output VCF (full path)')

# parser.add_argument('integers', metavar='N', type=int, nargs='+',
# 					help='an integer for the accumulator')

# parser.add_argument('--sum', dest='accumulate', action='store_const',
# 					const=sum, default=max,
# 					help='sum the integers (default: find the max)')

# print args.accumulate(args.integers)

args = parser.parse_args()

input_vcf = args.input_vcf[0]
output_vcf = os.path.dirname(input_vcf) + "/" + os.path.basename(input_vcf)[:-3] + "deleterious_only.vcf"

vcf_reader = vcf.Reader(open(input_vcf, 'r'))
vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)
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
	FILTER = ",".join(record.FILTER)


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



	if impact == 'HIGH' and LoF == "HC" and (ExAC_AF_NFE == "" or float(ExAC_AF_NFE) <= 0.01):
		pass_filter = True
	elif impact == 'MODERATE' and mut_type.find('missense') != -1 and (CADD_PHRED == "" or float(CADD_PHRED) > 30) and (ExAC_AF_NFE == "" or float(ExAC_AF_NFE) <= 0.01):
		pass_filter = True
	else:
		pass_filter = False

	# check also if variant is monomorphic
	n = len(record.samples)
	is_monomorphic = (record.num_het == n or record.num_hom_ref == n or record.num_hom_alt == n)
	pass_filter = pass_filter and not is_monomorphic and FILTER == ""

	if pass_filter:
		# print record.CHROM + '\t' + str(record.POS) + '\t' + str(record.ID) + '\t' + record.REF + '\t' + str(record.ALT[0]) + '\tFILTER=' + FILTER + '\tTYPE=' + mut_type + '\tIMPACT=' + impact + '\tGENE=' + gene + '\tCLIN_SIG=' + clin_sig + '\tLoF=' + LoF + '\tPosition=' + position + '\tExAC_MAF=' + ExAC_AF_NFE + '\tgnomAD_exomes_AF_nfe=' + gnomAD_exomes_AF_nfe + '\tCADD_PHRED=' + CADD_PHRED + '\tREVEL=' + REVEL + "\tHET=" + l_het + "\tHOM_ALT=" + l_hom_alts
		vcf_writer.write_record(record)



# python filter_variants.py /path_to_input_vcf/annotated.vcf



