#!/usr/bin/env python
import argparse
import sys
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser(description="")
parser.add_argument("-b", "--beds", nargs="+", help="bed file[s]")
parser.add_argument("-o", "--out", help="bed file representing regions with 1to1 mappings")
parser.add_argument("-w", "--window", help="window to check for overlaps in [100]", type=int, default=100)
args = parser.parse_args()

contigs = set()

def bed_to_tree(bedfile):
	tree = {}
	for line in open(bedfile):
		line = line.strip().split()
		contig, start, end = line[0], int(line[1]), int(line[2])
		if(contig not in tree):
			contigs.add(contig)
			tree[contig] = IntervalTree()
		tree[contig].add( Interval(start, end) )
	return(tree)

def read_genome(genome_file):
	genome = {}
	for line in open(genome_file):
		line = line.strip().split()
		contig, length = line[0], int(line[1])
		genome[contig] = length 
	return(genome)

# regions to test for 1 to 1
def to_test(trees):
	union = {}
	total = 0
	for contig in contigs:
		# all trees have this contig name
		if( sum( [contig in tree for tree in trees] ) == len(trees) ):
			union[contig] =  trees[0][contig].copy()
			for tree in trees:
				union[contig] |= tree[contig]
			union[contig].merge_overlaps()
			for interval in union[contig]:
				total += interval[1] - interval[0]
	return(union, total/args.window)



def is_1to1(contig, start, end, trees):
	for tree in trees:
		if(contig not in tree):
			return(False)
		overlaps = tree[contig][start:end]
		if(len(overlaps) != 1):
			return(False)
		# now only one overlap for sure
		o = overlaps.pop()
		o_st = o[0]
		o_en = o[1]
		# return false if overlap is not complete 
		if( o_st > start or o_en < end  ):
			return(False)
	return(True)


trees = [ bed_to_tree(bed) for bed in args.beds ]
valid, total = to_test(trees)

counter = 0
for contig in valid:
	intervals = [ (interval[0], interval[1]) for interval in valid[contig] ]
	for start, end in intervals:
		for cur_start in range(start, end, args.window):
			cur_end= min(end, cur_start + args.window) 
			if(not is_1to1(contig, cur_start, cur_end, trees)):
				valid[contig].chop(cur_start, cur_end)
			counter += 1
		sys.stderr.write("\rProgress:{:.3%}".format(counter/total))
sys.stderr.write("\n")


out = open(args.out, "w+")
for contig in sorted(valid):
	for inter in sorted(valid[contig]):
		out.write("{}\t{}\t{}\n".format(contig, inter[0], inter[1]))
out.close()

