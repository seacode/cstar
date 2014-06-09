# Path to ADMB distribution directory
export ADMB_HOME=/Users/stevenmartell1/admb-trunk/build/dist

# C compiler
CXX = clang++ 
CXXFLAGS = -O3 -DOPT_LIB -g -D__GNUDOS__ -Dlinux -DUSE_LAPLACE  -I. -I$(ADMB_HOME)/include -I$(ADMB_HOME)/contrib/include
LDFLAGS  = -shared  $(ADMB_HOME)/lib/libadmbo.a $(ADMB_HOME)/contrib/lib/libcontribo.a # linking flags
RM = rm -f  # rm command
TARGET_LIB = cstar.so # target lib
 
SRCS = dgamma.cpp dmultifan.cpp dnbinom.cpp dpois.cpp \
       generic.cpp 
       #approx.cpp selex.cpp
OBJS = $(SRCS:.cpp=.o)
 
.PHONY: all
all: ${TARGET_LIB}
 
$(TARGET_LIB): $(OBJS)
	@$(CXX) ${LDFLAGS} -o $@ $^
 
$(SRCS:.cpp=.d):%.d:%.cpp
	@echo 'compiling ' $<
	@$(CXX) $(CXXFLAGS) -MM $< >$@
 
include $(SRCS:.cpp=.d)
 
.PHONY: clean
clean:
	-${RM} ${TARGET_LIB} ${OBJS} $(SRCS:.cpp=.d)