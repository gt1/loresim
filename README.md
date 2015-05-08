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
