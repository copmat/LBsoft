#define CYCLE_OUT_INTERVAL(x, a,b) if ((x)<(a) .or. (x)>(b)) cycle

#define noALLAFAB
#define DEOWERN

#define noDEBUG_FORCEINT
#define noQUAD_FORCEINT
#define noLOGFORCES


#define IPRC 4
#define LPRC 4
#define PRC 8

#define D3Q19

#define noONLYCOM

#define NDIAGNSTREAM 1

#define NDIAGNHVAR 100

#define SVANBERG

#define IOOUT 6

#define _LB_(content) \
    CALL start_timing2("LB",#content);\
    CALL content;\
    CALL end_timing2("LB",#content)


#define RED(x) achar(27)//'[91m'//x//achar(27)//'[0m'
#define GREEN(x) achar(27)//'[92m'//x//achar(27)//'[0m'
#define YELLOW(x) achar(27)//'[93m'//x//achar(27)//'[0m'
#define PINK(x) achar(27)//'[95m'//x//achar(27)//'[0m'

#if PRC==4

#define C_MYFLOAT c_float

#define ZERO      0.e0
#define ONE       1.e0
#define TWO       2.e0
#define THREE     3.e0
#define FOUR      4.e0
#define FIVE      5.e0
#define SIX       6.e0
#define EIGHT     8.e0
#define NINE      9.e0
#define TEN      10.e0
#define ELEVEN   11.e0
#define TWELVE   12.e0
#define SIXTEEN  16.e0
#define EIGHTEEN  18.e0
#define TWENTYFOUR 24.e0
#define TWENTYSEVEN 27.e0
#define THIRTY   30.e0
#define THIRTYSIX 36.e0
#define FIFTY     50.e0
#define HALF      0.5e0
#define FOURTH    0.25e0

#define MINDENS   1.e-8

#elif PRC==8

#define C_MYFLOAT c_double

#define ZERO      0.d0
#define ONE       1.d0
#define TWO       2.d0
#define THREE     3.d0
#define FOUR      4.d0
#define FIVE      5.d0
#define SIX       6.d0
#define EIGHT     8.d0
#define NINE      9.d0
#define TEN      10.d0
#define ELEVEN   11.d0
#define TWELVE   12.d0
#define SIXTEEN  16.d0
#define EIGHTEEN  18.d0
#define TWENTYFOUR 24.d0
#define TWENTYSEVEN 27.d0
#define THIRTY   30.d0
#define THIRTYSIX 36.d0
#define FIFTY     50.d0
#define HALF      0.5d0
#define FOURTH    0.25d0

#define MINDENS   1.d-8

#else

//#error "ERROR in specifying PRC"

#endif

#ifdef D3Q19

#define LATTICE 319

#else  

//#error "ERROR in specifying type of lattice"

#endif
