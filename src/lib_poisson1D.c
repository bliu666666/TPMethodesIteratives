/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    int ii, jj, kk;
    for (jj=0;jj<(*la);jj++){
        kk = jj*(*lab);
        if (*kv>=0){
            for (ii=0;ii< *kv;ii++){
                AB[kk+ii]=0.0;
            }
        }
        AB[kk+ *kv]=-1.0;
        AB[kk+ *kv+1]=2.0;
        AB[kk+ *kv+2]=-1.0;
    }
    AB[0]=0.0;
    if (*kv == 1) {AB[1]=0;}

    AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
    int ii, jj, kk;
    for (jj=0;jj<(*la);jj++){
        kk = jj*(*lab);
        if (*kv>=0){
            for (ii=0;ii< *kv;ii++){
                AB[kk+ii]=0.0;
            }
        }
        AB[kk+ *kv]=0.0;
        AB[kk+ *kv+1]=1.0;
        AB[kk+ *kv+2]=0.0;
    }
    AB[1]=0.0;
    AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
    int jj;
    RHS[0]= *BC0;
    RHS[(*la)-1]= *BC1;
    for (jj=1;jj<(*la)-1;jj++){
        RHS[jj]=0.0;
    }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
    int jj;
    double h, DELTA_T;
    DELTA_T=(*BC1)-(*BC0);
    for (jj=0;jj<(*la);jj++){
        EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
    }
}  

void set_grid_points_1D(double* x, int* la){
    int jj;
    double h;
    h=1.0/(1.0*((*la)+1));
    for (jj=0;jj<(*la);jj++){
        x[jj]=(jj+1)*h;
    }
}

double relative_forward_error(double* x, double* y, int* la) {
    double ndiff = 0.0, ny = 0.0;
    for (int i = 0; i < *la; i++) {
        double d = x[i] - y[i];
        ndiff += d * d;
        ny += y[i] * y[i];
    }
    return sqrt(ndiff) / sqrt(ny);
}

int indexABCol(int i, int j, int *lab){
    return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
    int i;
    *info = 0;  // Initialiser info à 0, indiquant aucune erreur

    // Parcourir chaque colonne de la matrice et effectuer la décomposition LU de la matrice tridiagonale
    for (i = 0; i < *n - 1; i++) {
        if (AB[*ku * (*lab) + i] == 0.0) {  // Vérifier si la diagonale principale est 0
            *info = i + 1;  // Si un pivot zéro se produit, renvoie la position d'erreur
            return *info;
        }

        // Calculer L et mettre à jour U
        AB[(*ku - 1) * (*lab) + i + 1] /= AB[*ku * (*lab) + i];  // L = A[i+1][i] / A[i][i]
        AB[*ku * (*lab) + i + 1] -= AB[(*ku - 1) * (*lab) + i + 1] * AB[(*ku + 1) * (*lab) + i];  // U = A[i+1][i+1] - L * A[i][i+1]
    }

    // Définir la valeur du dernier pivot
    if (AB[*ku * (*lab) + *n - 1] == 0.0) {
        *info = *n;  // Le pivot de la dernière ligne est 0
    }

    return *info;  // Renvoie des informations, 0 indique le succès
}

void print_GB_operator_colMajor_poisson1D(double* AB, int ldab, int la, const char* name)
{
    printf("\n=== %s (ldab=%d, la=%d) in col-major band ===\n", name, ldab, la);
    for(int j=0; j<la; j++){
        printf("col %2d:  ", j);
        for(int r=0; r<ldab; r++){
            printf("%6.2f ", AB[r + j*ldab]);
        }
        printf("\n");
    }
    printf("=== End of %s ===\n\n", name);
}