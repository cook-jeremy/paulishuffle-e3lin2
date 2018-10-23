/**
#define width 4
#define height 4

#define depth 5
#define gatesPerLayer 4

#define circuitSize (depth*gatesPerLayer+1)*sizeof(Gate)

typedef enum GateTypes {
    HADAMARD, // not implemented
    PHASE,
    ROOTX,
    ROOTY,
    TGATE,
    NOISE, // not implemented
    CPHASEDOWN,
    CPHASERIGHT,
    OUTSTRING // must be the last gate
} GateType;


typedef struct Gates {
    GateType type; 
    uint64_t data; // 1 wherever gate is to be applied
} Gate;

#define tallyArraySize 256        // Must be larger than twice the nr of T gates
**/
#define tallyType int   // bit string of length 64
