ACLOCAL_AMFLAGS=-I m4

noinst_HEADERS=

MANPAGES = 

man_MANS = ${MANPAGES}

EXTRA_DIST = ${MANPAGES}

bin_PROGRAMS = loresim genseq fastareformat

loresim_SOURCES = loresim.cpp
loresim_LDADD = ${LIBMAUS2LIBS}
loresim_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
loresim_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

genseq_SOURCES = genseq.cpp
genseq_LDADD = ${LIBMAUS2LIBS}
genseq_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
genseq_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

fastareformat_SOURCES = fastareformat.cpp
fastareformat_LDADD = ${LIBMAUS2LIBS}
fastareformat_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
fastareformat_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

