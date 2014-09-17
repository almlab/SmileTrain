#!/usr/bin/env python

'''
author: scott w olesen (swo@mit.edu)
'''

from Bio import SeqIO, Seq
import argparse, sys

if __name__ == '__main__':
	p = argparse.ArgumentParser(description='reverse complement a single sequence, fasta, or list of sequences')
	g = p.add_mutually_exclusive_group(required=True)
	g.add_argument('-f', '--fasta', help='input fasta filename')
	g.add_argument('-s', '--seq', help='input sequence')
	g.add_argument('-l', '--list', help='input list')
	p.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help='output of same format')
	args = p.parse_args()

	if args.fasta:
		for record in SeqIO.parse(args.fasta, 'fasta'):
			record.seq = record.seq.reverse_complement()
			SeqIO.write(record, args.output, 'fasta')
	elif args.seq:
		seq = Seq.Seq(args.seq).reverse_complement()
		args.output.write("{}\n".format(seq))
	elif args.list:
		with open(args.list) as f:
			for line in f:
				seq = Seq.Seq(line.strip()).reverse_complement()
				args.output.write("{}\n".format(seq))