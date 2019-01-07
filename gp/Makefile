CFLAGS:= -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -fPIC -I\"/usr/local/include\"
LDFLAGS:= -bundle -undefined dynamic_lookup -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -fPIC


all: hyperellperiods.so

test: all
	gp < run.gp

%.so: %.o
	$(CC) -o $@ ${LDFLAGS} $<

#compile-command: "/usr/bin/gcc -c -o periods.gp.o -O3 -Wall
# * -fno-strict-aliasing -fomit-frame-pointer -fPIC -I\"/usr/local/include\"
#  *  periods.gp.c && /usr/bin/gcc -o periods.gp.so -bundle -undefined
#   *  dynamic_lookup -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -fPIC
#    *  periods.gp.o


