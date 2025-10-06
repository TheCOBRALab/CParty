#ifndef CONSTANTS_H
#define CONSTANTS_H

#define NONE 'N'  // no structure
#define HAIRP 'H' // closes a hairpin loop
#define INTER 'I' // closes an internal loop
#define MULTI 'M' // closes a regular multi-loop

#define M_WM 'B'  // closes a regular partial multi-loop
#define M_WMv 'v' // closes a regular partial multi-loop
#define M_WMp 'p' // closes a regular partial multi-loop

#define FREE 'W' // this base is free to be paired or unpaired to other base
#define LOOP 'V' // closes a loop

#define P_WMB 'R'
#define P_VP 'D'
#define P_WI 'G'
#define P_BE 'J'
#define P_WIP 'L'
#define P_WMBP 'T'
#define P_V 'A' // This case is only for the cases that we have some pairings in input structure that are not valid in terms of simfold restrictions
#define P_WMBW 'X'
#define P_VPL 'Y'
#define P_VPR 'Z'

#endif
