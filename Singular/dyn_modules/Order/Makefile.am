ACLOCAL_AMFLAGS = -I ../../m4

SOURCES = nforder.cpp nforder.h \
	  singular.cc \
	  nforder_elt.cc nforder_elt.h \
	  nforder_ideal.cc nforder_ideal.h 

MYINCLUDES =  -I${top_srcdir} -I${top_builddir} -I${srcdir} \
-I${top_srcdir}/libpolys -I${top_builddir}/libpolys \
-I${top_srcdir}/factory/include -I${top_builddir}/factory/include \
$(FACTORY_CFLAGS) $(NTL_CFLAGS) $(GMP_CFLAGS) $(GOOGLE_PERFTOOLS_CFLAGS)

module_LTLIBRARIES=nforder.la 

if ENABLE_P_PROCS_DYNAMIC
  moduledir = $(libexecdir)/singular/MOD
  P_PROCS_CPPFLAGS_COMMON = -DDYNAMIC_VERSION 
# Add under Mac OS X: -flat_namespace -weak_reference_mismatches weak -undefined dynamic_lookup
else
  moduledir = $(libdir)/singular
  P_PROCS_CPPFLAGS_COMMON = -DSTATIC_VERSION 
endif

P_PROCS_MODULE_LDFLAGS = -module -export-dynamic -avoid-version -flat_namespace -weak_reference_mismatches weak -undefined dynamic_lookup

nforder_la_SOURCES   = $(SOURCES)
nforder_la_CPPFLAGS  = ${MYINCLUDES} ${P_PROCS_CPPFLAGS_COMMON}
nforder_la_LDFLAGS   = ${AM_LDFLAGS} ${P_PROCS_MODULE_LDFLAGS} ${GOOGLE_PERFTOOL_LDFLAGS} 

# AM_COLOR_TESTS=always
#
# TESTS_ENVIRONMENT  = SINGULARPATH='${abs_top_builddir}/Singular/LIB:${abs_top_srcdir}/Singular/LIB:${abs_top_builddir}/libpolys/polys/.libs:${abs_top_builddir}/factory/gftables:${abs_builddir}/.libs:${abs_srcdir}'
# TESTS_ENVIRONMENT += SINGULAR_ROOT_DIR='${abs_top_builddir}'
# TESTS_ENVIRONMENT += SINGULAR_BIN_DIR='${abs_top_builddir}/Singular'
#
# TESTS=test_release.sh
#
# EXTRA_DIST = nforder.tst ederc.tst test_clear_enum.tst test.sh $(TESTS)
