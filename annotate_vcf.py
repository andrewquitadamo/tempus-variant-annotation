import requests
import json
import sys

vcf_file = sys.argv[1]

if len(sys.argv) <= 2:
    output_file = sys.stderr
else:
    output_file = open(sys.argv[2], 'w')

variant_consequences = {}
exac_allele_freq = {}

# Set the level of the consequences so they can be sorted.
# 0 - MODIFIER, 1 - LOW, 2 - MODERATE, 3 - HIGH
# Values based on https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
consequence_severity = {'3_prime_UTR_variant': 0,
                        '5_prime_UTR_variant': 0,
                        'intron_variant': 0,
                        'non_coding_transcript_exon_variant': 0,
                        'splice_region_variant': 1,
                        'synonymous_variant': 1,
                        'stop_retained_variant': 2,
                        'missense_variant': 2,
                        'initiator_codon_variant': 2,
                        'stop_lost': 3,
                        'stop_gained': 3,
                        'splice_donor_variant': 3,
                        'splice_acceptor_variant': 3}

# Loop through the VCF file to create a list that contains the variants in a format that the
# EXaC API accepts. The /rest/bulk/variant endpoint takes variants in the form of
# "[\"chr-position-ref-alt\",\"chr-position-ref-alt\"]"
with open(vcf_file, 'r') as vf:
    api_formatted_variants = []
    for line in vf:
        # Ignore header lines
        if line.startswith("#"):
            continue

        chrom, pos, _, ref, alt, *rest = line.strip().split('\t')

        # There is one variant that has a position of /. Obviously it can't be used to
        # search the API, so it is skipped
        if pos == "/":
            continue

        # A raw string is used so the \ are kept
        api_formatted_variants.append(r'\"' + "-".join((chrom, pos, ref, alt)) + r'\"')

url = 'http://exac.hms.harvard.edu/rest/bulk/variant'

# Loop through the variants in groups of 400. The EXaC API doesn't seem to have
# documentation on this, but 400 is just about the upper limit of variants that
# can be included in one bulk query. If too many are included no results are returned.
for i in range(0, len(api_formatted_variants), 400):
    data = r'"[' + ",".join(api_formatted_variants[i:i+400]) + r']"'
    # eval() is required so that the variants are formatted in the way the API expects
    r = requests.post(url, eval(data))
    exac_annotations = json.loads(r.text)

    for variant in exac_annotations:
        # Use 'chr-position' as the key for the exac_allel_freq and variant_consequences dictionaries
        position = "-".join(variant.split('-')[0:2])

        # Get the allele frequencies from EXaC
        if 'allele_freq' in exac_annotations[variant]['variant']:
            allele_freq = exac_annotations[variant]['variant']['allele_freq']
        else:
            allele_freq = 'NA'

        exac_allele_freq[position] = allele_freq

        # Get the variant consequence from EXaC
        if exac_annotations[variant]['consequence'] is None or exac_annotations[variant]['consequence'] == {}:
            variant_consequences[position] = 'NA'
        else:
            # Loop through consequences for the variants that have more than one consequence
            if len(exac_annotations[variant]['consequence']) > 1:
                most_severe_consequence_val = -1
                most_severe_consequence = None

                for consequence in exac_annotations[variant]['consequence']:
                    if consequence_severity[consequence] == most_severe_consequence_val:
                        # If there are multiple consequences of the same level keep the one
                        # that sorts first alphabetically
                        if consequence < most_severe_consequence:
                            most_severe_consequence = consequence

                    if consequence_severity[consequence] > most_severe_consequence_val:
                        most_severe_consequence_val = consequence_severity[consequence]
                        most_severe_consequence = consequence

                variant_consequences[position] = most_severe_consequence

            else:
                # Get the consequence for variants that only have one
                variant_consequences[position] = list(exac_annotations[variant]['consequence'].keys())[0]
                continue

# Loop through the VCF file again. Get the total read depth for each variant and
# the read depth for the alternate allele. Calculate the alternate allele read
# depth percentage. Add these to the INFO field along with the consequence and
# EXaC allele frequency.
with open(vcf_file, 'r') as vf:
    for line in vf:
        if line.startswith("#"):
            print(line, end='', file=output_file)

            if line.startswith("##INFO=<ID=END"):
                print('##INFO=<ID=ANNOT,Description="Variant annotations. Alternative allele read depth | Total read depth | Alternative allele read depth percentage | Variant consequence | EXaC allele frequency">', file=output_file)
            continue

        chrom, pos, ID, ref, alt, qual, Filter, info, Format, *genotypes = line.strip().split('\t')

        if pos == '/':
            continue
        # Separate the INFO field into its components. Extract the alternate allele depth (AO)
        # and total read depth (DP).
        info = info.split(';')
        alt_depth = info[5].split('=')[-1]
        read_depth = info[7].split('=')[-1]

        # Calculate the alternate allele read depth percentage. Perform the calculation
        # for each alternate allele in multiallelic variants.
        if ',' in alt_depth:
            depth_perc = []
            for depth in alt_depth.split(','):
                depth_perc.append('{:.3f}'.format((float(depth)/float(read_depth))))
            depth_perc = ','.join(depth_perc)
        else:
            depth_perc = '{:.3f}'.format(float(alt_depth)/float(read_depth))

        # Get the variants position to use as the key in variant_consequences and
        # exac_allele_freq dictionaries.
        position = chrom + "-" + pos

        # Add the annotation to the INFO field. Annotation fields are separated by |
        info = ';'.join(info) + ';ANNOT=' + ('|').join((alt_depth, read_depth, depth_perc, variant_consequences[position], str(exac_allele_freq[position])))
        print(chrom, pos, ID, ref, alt, qual, Filter, info, Format, '\t'.join(genotypes), sep='\t', file=output_file)
