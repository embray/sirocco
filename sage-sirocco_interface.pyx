
#### magic comment that tells Cython to link to libsirocco.so
#clib sirocco





cdef extern from "stdlib.h":
    void free(void* ptr)


# you can only include headers
cdef extern from "sirocco.h":
    double* homotopyPath (int degree, double *_coef, double _y0R, double _y0I)



def contpath(deg,values,y0r,y0i):
    cdef double* rop
    cdef double* c_values = <double*> sage_malloc(sizeof(double)*len(values))
    cdef int clen = <int> len(values)
    for i,v in enumerate(values):
        c_values[i] = values[i]
        
        
    
    cdef double y0R = y0r
    cdef double y0I = y0i
    sig_on()
    rop = homotopyPath (int(deg), c_values, y0R, y0I)
    sig_off()
    if rop == NULL:
        raise ValueError("libsirocco could not guarantee one step")
    n=int(rop[0])
    l=[0 for i in range(n)]
    for i in range(n):
        l[i]=(rop[3*i+1],rop[3*i+2],rop[3*i+3])
    free(rop)
    free(c_values)
    return l
#    return clen