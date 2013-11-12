CCTOOLS_HOME = /afs/crc.nd.edu/user/d/dpandiar/Public/cctools

PROGRAMS = elastic_sort 

LOCAL_LDFLAGS= -lwork_queue -ldttools -lm

all: ${PROGRAMS} 

elastic_sort: elastic_sort.c 
	gcc $^ -o $@ -I${CCTOOLS_HOME}/include/cctools -L${CCTOOLS_HOME}/lib ${LOCAL_LDFLAGS} 

clean:
	rm -f *~ *.o $(PROGRAMS)
