# tempus-variant-annotation

`annotate_vcf.py` takes a VCF file and annotates the variants with: Read depth at variant site, the number of reads with the alternate allele, the percentage of reads with the alternate allele, the consequence from EXaC and the allele frequency from EXaC. 

The input VCF file is supplied through a command line argument. A second argument specifies the output file to write the annotated VCF file to.
```shell
python3 annotate_vcf.py input.vcf output.vcf
```

If an output file is not specified the output is written to `STDOUT`, and can then be piped to any downstream scripts.
```shell
python3 annotate_vcf.py input.vcf
```

These annotations are appended to the `INFO` field using the `ANNOT` tag. The annotations are separated by a `|`. This format is based off of the SnpEff annotations http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf. 

`annotate_vcf.py` loops through the VCF file twice, once to extract the variant IDs to use in the EXaC API call and once to extract the read depth information and then print out the annotations.

The http://exac.hms.harvard.edu/rest/bulk/variant endpoint is used and groups of 400 variants are queried at a time. 400 seems to be close to the upper limit of what the API can handle, but the exact number does not seem to be documented. 

Currently only about half of the variants are annotated with consequence and EXaC allele frequency as they are not in the EXaC database. Another source for consequence predictions could be included in the future.
If there are multiple consequences found the most severe consequence is picked. If there are multiple consequences of the same severity than the first one alphabetically is picked. 

*Note:* The original VCF file included information on total read depth, and alternate allele read depth under the `DP` and `AO` tags in the `INFO` field. They have now also been included under the `ANNOT` tag as well. 

Dependancies:  
-------------- 
* Python3  
* requests - see http://docs.python-requests.org/en/master/user/install/#install for installation instructions  

Potential TODOs:
--------------------
[ ] Automatically populate the consequence_severity dictionary from an online table to include all possible consequences.  
[ ] Include a second annotation source like SnpEff or VEP  
[ ] Update the behavior of picking consequences with matching severity (Include all?)  
[ ] Could refactor into functions so unittests could be added  
