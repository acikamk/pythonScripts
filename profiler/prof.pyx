# prof.pyx
cdef extern from "/home/H0001-1/a.kovachev/pybios/gperftools/include/gperftools/profiler.h":
    void ProfilerStart( char* fname )
    void ProfilerStop()
 
def profiler_start(fname):
    ProfilerStart(<char *>fname)
 
def profiler_stop():
    ProfilerStop()