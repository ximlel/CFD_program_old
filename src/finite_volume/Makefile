CC = gcc
CFLAGS += -c
INC = ../include
AR = ar
ARFLAGS = crv
DIR_NAME = $(shell basename `pwd`)
FILES_C = $(shell ls *.c)

define compile_c_file
@for file in $(FILES_C); do \
( echo "$(CC) $(CFLAGS) $$file -I $(INC)" && $(CC) $(CFLAGS) $$file -I $(INC) ) \
done;
endef

define rm_o_file
@for file in $(FILES_O); do \
( rm -vf $$file ) \
done;
endef

all : compile archive

.PHONYP:all

compile:

	$(call compile_c_file)

.PHONYP:compile

FILES_O = $(shell ls *.o)

archive:

	$(AR) $(ARFLAGS)  $(DIR_NAME).a $(FILES_O)

./PHONYP:archive

clean:

	$(call rm_o_file)

.PHONYP:clean
