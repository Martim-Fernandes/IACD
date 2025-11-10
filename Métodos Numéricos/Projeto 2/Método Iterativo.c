#include <stdio.h>  
#include <math.h>   

#define EPSILON 1e-12 
#define MAX 1000 

double g1(double x) {
    return (x - pow(x, 3) - 4 * pow(x, 2) + 10);
}

double g2(double x) {
    return sqrt((10 / x) - 4 * x);  
}

double g3(double x) {
    return sqrt(10 - pow(x, 3)) * 1 / 2;  
}

double g4(double x) {
    return sqrt(10 / (4 + x));  
}

double g5(double x) {
    return (2 * pow(x, 3) + 4 * pow(x, 2) + 10) / (3 * pow(x, 2) + 8 * x);
}
void iterativo(double (*funcao)(double)) {
    double x0 = 1.5;  
    double x1;      
    int iter = 0;   
    do {
        x1 = funcao(x0); 
        iter++;
        printf("x%d = %.10f\n", iter, x1);
        if (iter >= MAX) {
            printf("Não convergiu em %d iterações\n", MAX);
            break;  
        }
        if (isnan(x1) || isinf(x1)) {
            printf("Divergiu\n");
            break;
        }
        if (fabs(x1 - x0) < EPSILON) {
            printf("Convergiu em %d iterações\n", iter);
            break;
        }
        x0 = x1;
    } while (1);  
    printf("\n");  
}

int main() {
    iterativo(g1);
    iterativo(g2);
    iterativo(g3);
    iterativo(g4);
    iterativo(g5);
    return 0;  
}
