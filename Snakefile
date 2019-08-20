import pandas as pd



DIR = os.path.dirname(workflow.snakefile) 
CWD = os.getcwd()



configfile:"asm.urls.txt"

sms = list(config.keys())

urls = {}
logs = {}
for sm in sms:
	urls[sm] = config[sm]["asm"]
	logs[sm] = config[sm]["log"]
	
asmfmt = "nobackups/asm/{sm}/{ID}.fasta"
bamfmt = "nobackups/asm/{sm}/{ID}.bam"
bedfmt = "nobackups/asm/{sm}/{ID}.bed"
conbed = "results/concensous.{sm}.bed",
quastbedfmt = "nobackups/asm/{sm}/{ID}.quast.bed"
filtquast = "results/{sm}.{ID}.quast.filt.bed"

wildcard_constraints:
	ID = "\d+",
	sm = "|".join(sms),

def get_final(notused):
	rtn = []
	for sm in sms:
		for ID in range(len( urls[sm] )):
			rtn.append( filtquast.format(sm=sm, ID=ID) )
	return(rtn)

def get_quast(notused):
	rtn = []
	for sm in sms:
		for ID in range(len( urls[sm] )):
			rtn.append( quastbedfmt.format(sm=sm, ID=ID) )
	return(rtn)


rule all:
	input:
		final = get_final,
		conbeds = expand(conbed, sm=sms),
		tbl = "results/summary.txt",

rule download:
	input:
	output:
		asm = asmfmt,
	resources:
		mem=4
	threads: 1
	run:
		url = urls[ wildcards["sm"] ][ int(str(wildcards["ID"])) ]	
		shell("wget -O {output.asm} {url}")

rule map:
	input:
		ref = "/net/eichler/vol26/projects/sda_assemblies/nobackups/assemblies/hg38/ucsc.hg38.no_alts.fasta", 
		asm = rules.download.output.asm,
	output:
		bam = bamfmt,
		bai = bamfmt + ".bai",
	resources:
		mem=4
	threads: 16
	shell:"""
minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5  {input.ref} {input.asm} | samtools view -F 260 -u - | samtools sort -m {resources.mem}G -@ {threads} - > {output.bam}  
samtools index {output.bam}
"""


rule bed:
	input:
		bam = rules.map.output.bam,
		bai = rules.map.output.bai,
	output:
		bed = bedfmt,
	resources:
		mem=16
	threads: 1
	shell:"""
bedtools bamtobed -i {input.bam} > {output.bed}
"""

def get_beds(wildcards):
	sm = str(wildcards.sm)
	rtn = []
	for ID in range(len( urls[sm] )):
		rtn.append( bedfmt.format(sm=sm, ID=ID) )
	return(rtn)

rule make_1_to_1_bed:
	input:
		beds = get_beds,
	output:
		bed = conbed, 
	resources:
		mem=16
	threads: 1
	shell:"""
{DIR}/1to1regions_with_windows.py -w 100 -b {input.beds}  -o {output.bed}
"""
	


def get_log(wildcards):
	sm = str(wildcards.sm)
	ID = int(str(wildcards.ID))
	#print(sm, ID, type(ID))
	#print(logs[sm])
	return(logs[sm][ID])


rule inter_quast:
	input:
		conbed = rules.make_1_to_1_bed.output.bed,
		quast = get_log,
	output:
		quast = quastbedfmt, 
		filtquast = filtquast, 
	resources:
		mem=16
	threads: 1
	shell:'''
grep "^chr" {input.quast} | grep -v "chrY" | awk '{{print $1"\t"$2"\t"$3"\t"$4 }}' | bedtools sort -i - > {output.quast}
bedtools intersect -f 1.0 -a {output.quast}  -b {input.conbed} > {output.filtquast} 
'''



rule summary:
	input:
		filt = get_final,
		bed = get_quast,
	output:
		tbl = "results/summary.txt",
	resources:
		mem=16
	threads: 1
	run:
		out = "Sample\tTotal\tConcensous\tFile\n"
		for filt, bed in zip(input.filt, input.bed):
			match = re.match("results/(.+).(\d+).quast.filt.bed", filt)
			Sample, ID = match.groups()
			ID=int(ID)
			Concensous = len( pd.read_csv(filt, sep="\t", header=None).index )
			Total = len( pd.read_csv(bed, sep="\t", header=None).index )
			File = config[Sample]["asm"][ID]
			
			out += f"{Sample}\t{Total}\t{Concensous}\t{File}\n"

		print(out)	
		open(output["tbl"], "w+").write(out)











