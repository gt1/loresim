loresim
=======

Long read simulation.

Source
------

The loresim source code is hosted on github:

	git@github.com:gt1/loresim.git

Compilation of loresim
----------------------

loresim needs libmaus2 [https://github.com/gt1/libmaus2] . When libmaus2
is installed in ${LIBMAUSPREFIX} then biobambam2 can be compiled and
installed in ${HOME}/biobambam2 using

	- autoreconf -i -f
	- ./configure --with-libmaus2=${LIBMAUSPREFIX} \
		--prefix=${HOME}/biobambam2
	- make install

Contained programs
------------------

*genseq*: generates a sequence of random DNA blocks in FastA format. Blocks are repeated with a given probability. The program has the following options

- blocksizeavg: average block size (default: 8k)
- blocksizestd: block size standard deviation (default: 4k)
- blockmin: minimum block size (default: 512)
- repfrag: block repeat probability (default: 0.01)
- reprepfrag: probability that repeated block is chosen from list of already repeated blocks (default: 0.5)
- numblocks: number of blocks to be produced (default: 1k)
- erate: rate of error spiked into output blocks (default: 0.001)
- eratesq: error rate standard deviation (default: 0.001)
- substrate: fraction of substitutions in introduced base errors
- delrate: fraction of deletions in introduced base errors
- insrate: fraction of insertions in introduced base errors

*fastareformat*: reformats FastA or FastQ so it becomes valid input for the DAZZ_DB tool fasta2DB. This includes changing read names and wrapping long lines.

*loresim*: generates a sequence of random DNA reads from a given sequence. For this the input sequence is traversed a number of times and split into pieces which are subsequently either
modified within a certain error range or dropped with a certain rate. The program can switch between a low and high error profile during read generation. Reads are randomly taken from the forward or reverse strand.
The program has the following options

- readlenavg: average read length (default: 15k)
- readlenstddev: standard deviation for read length (default: 3k)
- droprate: probability for dropping a read (default: 0.01)
- numtraversals: number of traversals of the input sequence (default: 1)
- startlowprob: probability of starting a read in low error mode (default: 0.7)
- keeplowstate: probability of staying in low error mode when in low error mode (default: 0.9998)
- keephighstate: probability of staying in high error mode when in high error mode (default: 0.995)
- eratelow: low error mode error rate average
- eratelowstddev: standard deviation of error rate in low error mode
- eratehigh: high error mode error rate average
- eratehighstddev: standard deviation of error rate in high error mode
