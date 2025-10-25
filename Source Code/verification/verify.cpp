#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[])
{
    char *name1 = argv[1];
    char *name2 = argv[2];

    FILE *file1 = fopen(name1, "r");
    FILE *file2 = fopen(name2, "r");

    double value1, value2;
    while (fscanf(file1, "%lf", &value1) == 1 && fscanf(file2, "%lf", &value2) == 1)
    {
        if (fabs(value1 - value2) > 1e-6)
        {
            printf("Values differ: %f vs %f\n", value1, value2);
        }        
    }
    printf("Comparison complete.\n");
    return 0;
}