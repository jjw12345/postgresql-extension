# contrib/recommend/Makefile

MODULE_big = recommend
OBJS =recommend.o smlar_gist.o smlar_cache.o smlar_guc.o \
       smlar_stat.o 
PGFILEDESC = "recommend - similarity data compute for PostgreSQL"


EXTENSION = recommend
DATA = recommend--1.0.sql

ifdef USE_PGXS

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
else
subdir = contrib/recommend
top_builddir = ../..
include $(top_builddir)/src/Makefile.global
include $(top_srcdir)/contrib/contrib-global.mk
endif
