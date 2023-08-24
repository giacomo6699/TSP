#include <stdio.h>
#include <stdlib.h>

int main(){
    FILE *out = fopen("./plot/data.dat", "w");
    fprintf(out, "0 0\r\n");
    fprintf(out, "1 1\r\n");
    fprintf(out, "2 0\r\n");
    fprintf(out, "3 1\r\n");
    fclose(out);
    system("gnuplot ./plot/commands.txt");
}