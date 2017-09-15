/**
 * Document: MaxCompiler Tutorial (maxcompiler-tutorial.pdf)
 * Chapter: 6      Example: 1      Name: Typecast
 * MaxFile name: TypeCast
 * Summary:
 *        Casts from an unsigned int to a float and back.
 */
#include <stdio.h>
#include <stdlib.h>
#include "Maxfiles.h"
#include <MaxSLiCInterface.h>
#include <time.h>

/*
 * This implementation will formulate the consensus string as the final output
 */
 
/*  Size of a DNA strand*/
#define SIZE 1048576
#define motifSIZE 4
#define A 65
#define T 84
#define G 71
#define C 67

/* basic 5 dna sequences */
static int8_t dna_seq1[SIZE],dna_seq2[SIZE],dna_seq3[SIZE],dna_seq4[SIZE],dna_seq5[SIZE];
/*  for storing profile matrix*/
static int8_t profile_A[SIZE],profile_T[SIZE],profile_G[SIZE],profile_C[SIZE];
/*  number of As,Ts,Gs and Cs*/

int nA = 0;
int nT = 0;
int nG = 0;
int nC = 0;
/*  for profile matrix in CPU mode(without dataflow)*/
int profile_CPU[4][SIZE];
/*  for counting stop codons found*/
static int stopCodons = 0;

int i,j,maxscore_CPU,maxscore_DFE;
int colMax,loc,score,temp,val;

/*
 * Loops in initialize() contribute significantly to the overall execution time
 *  a single loop is used again and again for various dna_sequences
 *  this code is just for a random DNA generator, so this won't be used in actual experiments.
 */
void initialize(){
    int val;
    time_t t;

/*  initialize the random number generator  */
    srand((unsigned) time(&t));

    for (i = 0; i < SIZE; i++) {
        val = rand() % 4;
        if (val == 0) dna_seq1[i] = A;
        else if (val == 1) dna_seq1[i] = T;
        else if (val == 2) dna_seq1[i] = G;
        else if (val == 3) dna_seq1[i] = C;
    }
    for (i = 0; i < SIZE; i++) {
        val = rand() % 4;
        if (val == 0) dna_seq2[i] = A;
        else if (val == 1) dna_seq2[i] = T;
        else if (val == 2) dna_seq2[i] = G;
        else if (val == 3) dna_seq2[i] = C;
    }
    for (i = 0; i < SIZE; i++) {
        val = rand() % 4;
        if (val == 0) dna_seq3[i] = A;
        else if (val == 1) dna_seq3[i] = T;
        else if (val == 2) dna_seq3[i] = G;
        else if (val == 3) dna_seq3[i] = C;
    }
    for (i = 0; i < SIZE; i++) {
        val = rand() % 4;
        if (val == 0) dna_seq4[i] = A;
        else if (val == 1) dna_seq4[i] = T;
        else if (val == 2) dna_seq4[i] = G;
        else if (val == 3) dna_seq4[i] = C;
    }
    for (i = 0; i < SIZE; i++) {
        val = rand() % 4;
        if (val == 0) dna_seq5[i] = A;
        else if (val == 1) dna_seq5[i] = T;
        else if (val == 2) dna_seq5[i] = G;
        else if (val == 3) dna_seq5[i] = C;
    }
}

/*  for calculating the number and displaying the locations of the stop codons*/
void stopcodons(int type,int dna,int loc){
    stopCodons++;
/*
    if(type == 1){
        printf("ATC (amber) located in DNA sequence = %d at location = %d\n",dna,loc);
    }
    else if(type == 2){
        printf("ATT (ochre) located in DNA sequence = %d at location = %d\n",dna,loc);
    }
    else if(type == 3){
        printf("ACT (umber) located in DNA sequence = %d at location = %d\n",dna,loc);
    }
    else if(type == 4){
        printf("TAG (amber) located in DNA sequence = %d at location = %d\n",dna,loc);
    }
    else if(type == 5){
        printf("TAA (ochre) located in DNA sequence = %d at location = %d\n",dna,loc);
    }
    else if(type == 6){
        printf("TGA (umber) located in DNA sequence = %d at location = %d\n",dna,loc);
    }
*/
}

/*
 * There are 5 loops in calculateProfile_CPU() which can be remodeled to be one
 * it contributes significantly to the execution time
 * will execute for every entry in the DNA sequence
 * has to deal with huge data i.e. 5 DNA strands
 *
 * CAN BE CONSIDERED FOR PORTING
 */
void calculateProfile_CPU() {
    for (i = 0; i < SIZE; i++) {
        if(dna_seq1[i]==A) nA++;
        else if(dna_seq1[i]==T) nT++;
        else if(dna_seq1[i]==G) nG++;
        else if(dna_seq1[i]==C) nC++;

        if(dna_seq2[i]==A) nA++;
        else if(dna_seq2[i]==T) nT++;
        else if(dna_seq2[i]==G) nG++;
        else if(dna_seq2[i]==C) nC++;

        if(dna_seq3[i]==A) nA++;
        else if(dna_seq3[i]==T) nT++;
        else if(dna_seq3[i]==G) nG++;
        else if(dna_seq3[i]==C) nC++;

        if(dna_seq4[i]==A) nA++;
        else if(dna_seq4[i]==T) nT++;
        else if(dna_seq4[i]==G) nG++;
        else if(dna_seq4[i]==C) nC++;

        if(dna_seq5[i]==A) nA++;
        else if(dna_seq5[i]==T) nT++;
        else if(dna_seq5[i]==G) nG++;
        else if(dna_seq5[i]==C) nC++;


        profile_CPU[0][i] = nA; nA = 0;
        profile_CPU[1][i] = nT; nT = 0;
        profile_CPU[2][i] = nG; nG = 0;
        profile_CPU[3][i] = nC; nC = 0;
    }
}

/*
 * Single nested loop is there in maxScore_CPU() which performs critical work
 * takes significant amount of execution time
 * used only once
 * deals with huge data of size = SIZE*4
 */
int maxScore_CPU(){
    int temp;
    int max = 0;
    int cons[SIZE];
    for(j = 0; j < SIZE; j++){
        temp = 0;
        for(i = 0; i < 4; i++){
            if(temp < profile_CPU[i][j]){
                temp = profile_CPU[i][j];
                if(i == 0)      cons[j] = A;
                else if(i == 1) cons[j] = T;
                else if(i == 2) cons[j] = G;
                else if(i == 3) cons[j] = C;
            }
        }
        max = max + temp;
    }

    printf("Score for the given strands (without maxeler) : %d\n\n",max);
/*    printf("Consensus string : ");
    for(i = 0; i < SIZE; i++)
        printf("%c ",cons[i]);
    printf("\n\n");
*/
    return max;
}

int calcMax_CPU(){
    int max = 0;
    for(i = 0; i < SIZE - motifSIZE; i++){
        score = 0;
        for(j = 0; j < motifSIZE; j++) {
            temp = 0;
            for(int k = 0; k < 4; k++){
                if(temp < profile_CPU[k][i+j]){
                    temp = profile_CPU[k][i+j];
                }
            }
            score = score + temp;
        }
        if(score > max){
            max = score;
            loc = i;
        }
    }
    maxscore_CPU = max;
    return loc;
}

void consensus_CPU(int loc){
    int str[motifSIZE];
    int k = 0;
    for(i = loc; i < motifSIZE + loc; i++){
        colMax = 0;
        if(colMax < profile_CPU[0][i]){
            val = A;
            colMax = profile_CPU[0][i];
        }
        if(colMax < profile_CPU[1][i]){
            val = T;
            colMax = profile_CPU[1][i];
        }
        if(colMax < profile_CPU[2][i]){
            val = G;
            colMax = profile_CPU[2][i];
        }
        if(colMax < profile_CPU[3][i]){
            val = C;
            colMax = profile_CPU[3][i];
        }
        str[k++] = val;
    }

    printf("\nConsensus : ");
    for(i = 0; i < motifSIZE; i++) printf("%c ",str[i]);
    printf("\n");
}

void printDNAs(){
    printf("Given DNA strands are as follows :\n");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq1[i]);
    printf("\n");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq2[i]);
    printf("\n");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq3[i]);
    printf("\n");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq4[i]);
    printf("\n");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq5[i]);
    printf("\n\n");
}

void printProfile_CPU(){
    printf("Profile Matrix for CPUCode : \n");
    for(i = 0; i < 4; i++){
        for(j = 0; j < SIZE; j++) printf("%d ",profile_CPU[i][j]);
        printf("\n");
    }
    printf("\n");
}

void printProfile_DFE(){
    printf("Profile Matrix for DFECode : \n");
    for(i = 0; i < SIZE; i++) printf("%d ",profile_A[i]);
    printf("\n");
    for(i = 0; i < SIZE; i++) printf("%d ",profile_T[i]);
    printf("\n");
    for(i = 0; i < SIZE; i++) printf("%d ",profile_G[i]);
    printf("\n");
    for(i = 0; i < SIZE; i++) printf("%d ",profile_C[i]);
    printf("\n\n");
}

int CPU(){
    int loc;
    calculateProfile_CPU();
    printProfile_CPU();

    printf("\nCPU Method\n");
    loc = calcMax_CPU();
    printf("\nloc = %d\n",loc+1);
    printf("\nMax Score(s,DNA) = %d\n",maxscore_CPU);

    consensus_CPU(loc);
    return maxscore_CPU;
}
int calcMax_DFE(){
    maxscore_DFE = 0;
    for(i = 0; i < SIZE - motifSIZE; i++){
        score = 0;
        for(j = 0; j < motifSIZE; j++) {
            colMax = 0;
            if(colMax < profile_A[i+j]) colMax = profile_A[i+j];
            if(colMax < profile_T[i+j]) colMax = profile_T[i+j];
            if(colMax < profile_G[i+j]) colMax = profile_G[i+j];
            if(colMax < profile_C[i+j]) colMax = profile_C[i+j];
            score = score + colMax;
        }
        if(maxscore_DFE < score){
            maxscore_DFE = score;
            loc = i;
        }
    }
    return loc;
}

void consensus_DFE(int loc){
    int str[motifSIZE];
    int k = 0;
    for(i = loc; i < motifSIZE + loc; i++){
        colMax = 0;
        if(colMax < profile_A[i]){
            val = A;
            colMax = profile_A[i];
        }
        if(colMax < profile_T[i]){
            val = T;
            colMax = profile_T[i];
        }
        if(colMax < profile_G[i]){
            val = G;
            colMax = profile_G[i];
        }
        if(colMax < profile_C[i]){
            val = C;
            colMax = profile_C[i];
        }
        str[k++] = val;
    }

    printf("\nConsensus : ");
    for(i = 0; i < motifSIZE; i++) printf("%c ",str[i]);
    printf("\n");
}

/*
 *  this function will be totally to DFE
 */
int DFE(){
    
    TypeCast(SIZE,dna_seq1,dna_seq2,dna_seq3,dna_seq4,dna_seq5,profile_A,profile_T,profile_G,profile_C);

    printProfile_DFE();

    printf("\n\nDFE Method\n");

    loc = calcMax_DFE();
    printf("\nloc = %d\n",loc+1);
    printf("\nMax Score(s,DNA) = %d\n",maxscore_DFE);

    consensus_DFE(loc);

    return maxscore_DFE;
}

void checkStopCodons(){
    for(i = 0; i < SIZE; i++) {
        if (dna_seq1[i] == A) {
            if (dna_seq1[i + 1] == T) {
                if (dna_seq1[i + 2] == C) stopcodons(1, 1, i);
                else if (dna_seq1[i + 2] == T) stopcodons(2, 1, i);
            } else if (dna_seq1[i + 1] == C && dna_seq1[i + 2] == T) stopcodons(3, 1, i);
        } else if (dna_seq1[i] == T) {
            if (dna_seq1[i + 1] == A) {
                if (dna_seq1[i + 2] == G) stopcodons(4, 1, i);
                else if (dna_seq1[i + 2] == A) stopcodons(5, 1, i);
            } else if (dna_seq1[i + 1] == G && dna_seq1[i + 2] == A) stopcodons(6, 1, i);
        }

        if (dna_seq2[i] == A) {
            if (dna_seq2[i + 1] == T) {
                if (dna_seq2[i + 2] == C) stopcodons(1, 2, i);
                else if (dna_seq2[i + 2] == T) stopcodons(2, 2, i);
            } else if (dna_seq2[i + 1] == C && dna_seq2[i + 2] == T) stopcodons(3, 2, i);
        } else if (dna_seq2[i] == T) {
            if (dna_seq2[i + 1] == A) {
                if (dna_seq2[i + 2] == G) stopcodons(4, 2, i);
                else if (dna_seq2[i + 2] == A) stopcodons(5, 2, i);
            } else if (dna_seq2[i + 1] == G && dna_seq2[i + 2] == A) stopcodons(6, 2, i);
        }

        if (dna_seq3[i] == A) {
            if (dna_seq3[i + 1] == T) {
                if (dna_seq3[i + 2] == C) stopcodons(1, 3, i);
                else if (dna_seq3[i + 2] == T) stopcodons(2, 3, i);
            } else if (dna_seq3[i + 1] == C && dna_seq3[i + 2] == T) stopcodons(3, 3, i);
        } else if (dna_seq3[i] == T) {
            if (dna_seq3[i + 1] == A) {
                if (dna_seq3[i + 2] == G) stopcodons(4, 3, i);
                else if (dna_seq3[i + 2] == A) stopcodons(5, 3, i);
            } else if (dna_seq3[i + 1] == G && dna_seq3[i + 2] == A) stopcodons(6, 3, i);
        }

        if (dna_seq4[i] == A) {
            if (dna_seq4[i + 1] == T) {
                if (dna_seq4[i + 2] == C) stopcodons(1, 4, i);
                else if (dna_seq4[i + 2] == T) stopcodons(2, 4, i);
            } else if (dna_seq4[i + 1] == C && dna_seq4[i + 2] == T) stopcodons(3, 4, i);
        } else if (dna_seq4[i] == T) {
            if (dna_seq4[i + 1] == A) {
                if (dna_seq4[i + 2] == G) stopcodons(4, 4, i);
                else if (dna_seq4[i + 2] == A) stopcodons(5, 4, i);
            } else if (dna_seq4[i + 1] == G && dna_seq4[i + 2] == A) stopcodons(6, 4, i);
        }

        if (dna_seq5[i] == A) {
            if (dna_seq5[i + 1] == T) {
                if (dna_seq5[i + 2] == C) stopcodons(1, 5, i);
                else if (dna_seq5[i + 2] == T) stopcodons(2, 5, i);
            } else if (dna_seq5[i + 1] == C && dna_seq5[i + 2] == T) stopcodons(3, 5, i);
        } else if (dna_seq5[i] == T) {
            if (dna_seq5[i + 1] == A) {
                if (dna_seq5[i + 2] == G) stopcodons(4, 5, i);
                else if (dna_seq5[i + 2] == A) stopcodons(5, 5, i);
            } else if (dna_seq5[i + 1] == G && dna_seq5[i + 2] == A) stopcodons(6, 5, i);
        }
    }
    printf("\nSTOP codons found = %d\n\n\n",stopCodons);
}

int main(){

    printf("Simplest Method for solving Motif Finding Problem\n");
    printf("Implemented by Jayendra Khandare(jkhandar@umail.iu.edu)\t DATE : June 13,2017\n\n");

    int cpu,dfe;

    initialize();
    printDNAs();

    cpu = CPU();

    dfe = DFE();

    if(cpu == dfe)  printf("\nSTATUS : OK\n");
    
    checkStopCodons();
    
    return 0;
}
