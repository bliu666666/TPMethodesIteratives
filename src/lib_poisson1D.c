/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    int i;
    for (i = 0; i < *la; i++) {
        // diagonale de la bande inférieure
        if (i > 0) AB[(*kv - 1) * (*la) + i] = -1.0;
        // diagonale principale
        AB[*kv * (*la) + i] = 2.0;
        // Diagonale de la bande supérieure
        if (i < *la - 1) AB[(*kv + 1) * (*la) + i] = -1.0;
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
    int i;
    for (i = 0; i < *la; i++) {
        // La diagonale de la bande inférieure et la diagonale de la bande supérieure sont définies sur 0
        if (i > 0) AB[(*kv - 1) * (*la) + i] = 0.0;
        AB[*kv * (*la) + i] = 1.0; // Diagonale principale réglée à 1
        if (i < *la - 1) AB[(*kv + 1) * (*la) + i] = 0.0;
    }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
    int i;
    for (i = 0; i < *la; i++) {
        RHS[i] = 0.0;
    }
    RHS[0] += *BC0;
    RHS[*la - 1] += *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
    int i;
    for (i = 0; i < *la; i++) {
        EX_SOL[i] = (*BC0) + X[i] * ((*BC1) - (*BC0));
    }
}  

void set_grid_points_1D(double* x, int* la){
    int i;
    double h = 1.0 / ((double)(*la + 1)); // espace de la grille
    for (i = 0; i < *la; i++) {
        x[i] = (i + 1) * h; // Commencer par le premier point intérieur
    }
}

double relative_forward_error(double* x, double* y, int* la){
    double norm_diff = 0.0, norm_y = 0.0;
    int i;
    for (i = 0; i < *la; i++) {
        norm_diff += (x[i] - y[i]) * (x[i] - y[i]);
        norm_y += y[i] * y[i];
    }
    return sqrt(norm_diff) / sqrt(norm_y);
}

int indexABCol(int i, int j, int *lab){
    return 0;
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
