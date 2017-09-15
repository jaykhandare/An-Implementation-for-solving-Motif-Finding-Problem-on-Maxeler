#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 20
#define motifSIZE 5
#define A 65
#define T 84
#define G 71
#define C 67

static int8_t dna_seq1[SIZE],dna_seq2[SIZE],dna_seq3[SIZE],dna_seq4[SIZE],dna_seq5[SIZE];
static int8_t profile_A[SIZE],profile_T[SIZE],profile_G[SIZE],profile_C[SIZE];
static int8_t profile_CPU[4][SIZE];
int i,j,maxscore_CPU,maxscore_DFE;
int nA = 0;
int nT = 0;
int nG = 0;
int nC = 0;
int stopCodons = 0;

void initialize(){
    int val;
    time_t t;
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

void calculateProfile_DFE() {
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

        profile_A[i] = nA; nA = 0;
        profile_T[i] = nT; nT = 0;
        profile_G[i] = nG; nG = 0;
        profile_C[i] = nC; nC = 0;
    }
}

void printDNAs(){
    printf("\n5 DNA strands are generated as follows :\n");
    printf("D1 ");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq1[i]);
    printf("\n");
    printf("D2 ");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq2[i]);
    printf("\n");
    printf("D3 ");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq3[i]);
    printf("\n");
    printf("D4 ");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq4[i]);
    printf("\n");
    printf("D5 ");
    for(i = 0; i < SIZE; i++) printf("%c ",dna_seq5[i]);
}

void printProfile_DFE(){
    printf("\n\nProfile Matrix for DFE Code : \n");
    printf("A ");
    for(j = 0; j < SIZE; j++) printf("%d ",profile_A[j]);
    printf("\n");
    printf("T ");
    for(j = 0; j < SIZE; j++) printf("%d ",profile_T[j]);
    printf("\n");
    printf("G ");
    for(j = 0; j < SIZE; j++) printf("%d ",profile_G[j]);
    printf("\n");
    printf("C ");
    for(j = 0; j < SIZE; j++) printf("%d ",profile_C[j]);
    printf("\n");
}

int calcMax_DFE(){
    int colMax,loc,score;
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
    int colMax,val;
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

int DFE(){
    int loc;

    calculateProfile_DFE();
    printProfile_DFE();

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
    printf("\nSTOP codons found = %d\n",stopCodons);
}

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
void printProfile_CPU(){
    printf("\n\nProfile Matrix for CPUCode : \n");
    printf("A ");
    for(j = 0; j < SIZE; j++) printf("%d ",profile_CPU[0][j]);
    printf("\n");
    printf("T ");
    for(j = 0; j < SIZE; j++) printf("%d ",profile_CPU[1][j]);
    printf("\n");
    printf("G ");
    for(j = 0; j < SIZE; j++) printf("%d ",profile_CPU[2][j]);
    printf("\n");
    printf("C ");
    for(j = 0; j < SIZE; j++) printf("%d ",profile_CPU[3][j]);
    printf("\n");
}

void consensus_CPU(int loc){
    int str[motifSIZE];
    int colMax,val;
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

int calcMax_CPU(){
    int temp,loc,score;
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

int CPU(){
    int loc;
    calculateProfile_CPU();
    printProfile_CPU();
    loc = calcMax_CPU();
    printf("\nloc = %d\n",loc+1);
    printf("\nMax Score(s,DNA) = %d\n",maxscore_CPU);

    consensus_CPU(loc);

    return maxscore_CPU;
}

int main(){

    printf("Simplest Method for solving 'Motif Finding Problem'\n");
    printf("Implemented by Jayendra Khandare(jkhandar@umail.iu.edu)\t DATE : June 13,2017\n\n");

    initialize();
    printDNAs();

    printf("***CPU code***\n\n");
    int cpu = CPU();

    printf("***DFE code***\n\n");
    int dfe = DFE();

    checkStopCodons();

    if(cpu == dfe) printf("\nStatus : OK\n");

    return 0;
}
