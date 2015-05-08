#! /bin/bash
if [ $# -lt 3 ] ; then
	echo "usage: $0 sequence_in.fasta reads_out.fasta #traversals"
	exit 1
fi

if [ ! -f src/loresim ] ; then
	echo "Program src/loresim does not exist. Please compile it first."
	exit 1
fi

IN=$1
OUT=$2
TRAVERSALS=$3

OUTSHORT=${OUT%.fasta}
OUTSHORT=${OUTSHORT%.fa}
OUTSHORT=${OUTSHORT}_short_names.fasta

src/loresim numtraversals=$TRAVERSALS <$IN >$OUT
awk < $OUT '/^>/ { print $1 " " $2 } !/^>/ { print }' >${OUTSHORT}

if [ -f ../DAZZ_DB/fasta2DB ] ; then
	if [ -f ${OUTSHORT%.fasta}.db ] ; then
		../DAZZ_DB/DBrm ${OUTSHORT%.fasta}.db
	fi
	../DAZZ_DB/fasta2DB ${OUTSHORT%.fasta}.db ${OUTSHORT}

	if [ -f ../DALIGNER/HPCdaligner ] ; then
		rm -f ${OUTSHORT%.fasta}*.las
		../DALIGNER/HPCdaligner ${OUTSHORT%.fasta}.db | PATH=$PATH:../DALIGNER bash
		
		if [ -f ../../celamy2/LAcelamy ] ; then
			../../celamy2/LAcelamy ${OUTSHORT%.fasta}.db ${OUTSHORT%.fasta}.db out_short_names.1.las
		fi
	fi
fi
