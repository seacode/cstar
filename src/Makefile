# export ADMB_HOME=/Users/stevenmartell1/admb-trunk/build/dist
export ADMB_HOME="C:/Program Files (x86)/ADMB"

SRCS = dgamma.cpp dmultifan.cpp dnbinom.cpp dpois.cpp \
       generic.cpp
       #approx.cpp selex.cpp
OBJS = $(SRCS:.cpp=.obj)

all: clean
	$(MAKE) libcstar.a

libcstar.a: $(SRCS)
	$(ADMB_HOME)/admb -c $^
	$(AR) -rs  $@ $(OBJS)

clean:
	rm -vf ${TARGET_LIB} $(OBJS)
