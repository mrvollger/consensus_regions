# consensus_regions


## generate alignments against hg38 
```
minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5  {input.ref} {input.asm} | samtools view -F 260 -u - | samtools sort -m 8G -@ {threads} - > {output.bam}  

bedtools bamtobed -i {output.bam} | bedtools sort -i - | cut -f 1,2,3,4,5 > {output.bed}

```


## generate 1 to 1 regions
```
1to1regions_with_windows.py -w 100 --beds 1.bed 2.bed ... n.bed  -o consensus.bed 
```


## example comparision commands 
```

# input beds to filter
grep "^chr" ../CHM13_Shasta_all_missassemblies_total_1107.log | grep -v "chrY" | awk '{print $1"\t"$2"\t"$3"\t"$4 }' | bedtools sort -i - > shasta.chm13.log.bed
grep "^chr" ../CHM13_HiFi_all_missassemblies_total_8666.log | grep -v "chrY" | awk '{print $1"\t"$2"\t"$3"\t"$4 }' | bedtools sort -i - > canu.chm13.log.bed


# filter inputs by concenous regions in this case "consensus.bed" is  "hifi_vs_shasta_2_w100.bed"
bedtools intersect -f 1.0 -a shasta.chm13.log.bed  -b hifi_vs_shasta_2_w100.bed > shasta.chm13.1to1.log.bed
bedtools intersect -f 1.0 -a canu.chm13.log.bed  -b hifi_vs_shasta_2_w100.bed > canu.chm13.1to1.log.bed 

# count remaining entries 
wc -l  shasta.chm13.1to1.log.bed  canu.chm13.1to1.log.bed 

# intersect against eachother 
bedtools intersect -f 0.5 -r -a shasta.chm13.1to1.log.bed -b canu.chm13.1to1.log.bed  | wc -l 
bedtools intersect -f 0.9 -r -a shasta.chm13.1to1.log.bed -b canu.chm13.1to1.log.bed  | wc -l 

# intersect against know SVs
bedtools intersect -f 0.5 -r -a shasta.chm13.1to1.log.bed -b CHM13.SVs.bed  | wc -l 
bedtools intersect -f 0.5 -r -a canu.chm13.1to1.log.bed -b CHM13.SVs.bed  | wc -l 
```

