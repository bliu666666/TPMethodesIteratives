/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
    double h = 1.0 / ((*la) + 1);
    for (int i = 0; i < *la; i++) {
        eigval[i] = 4 * sin((i + 1) * M_PI * h / 2.0)*sin((i + 1) * M_PI * h / 2.0);
    }
}

double eigmax_poisson1D(int *la) {
    double h = 1.0 / ((*la) + 1);
    return 4 * sin((*la) * M_PI * h / 2.0)*sin((*la) * M_PI * h / 2.0);
}

double eigmin_poisson1D(int *la){
    double h = 1.0 / ((*la) + 1);
    return 4 * sin(M_PI * h / 2.0)*sin(M_PI * h / 2.0);
}

double richardson_alpha_opt(int *la){
    double lam_max = eigmax_poisson1D(la);
    double lam_min = eigmin_poisson1D(la);
    double alpha_opt = 2.0 / (lam_max + lam_min);
    return alpha_opt;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    int n = *la;
    double alpha = *alpha_rich;

    double *res = (double*) malloc(n*sizeof(double));
    int k=0;

    double resn=1.0;

    while((resn>*tol) && (k<(*maxit)-1)){
        // new res
        cblas_dcopy(n, RHS,1,res,1);
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
                    n, n, *kl, *ku,
                    -1.0, AB,*lab,
                    X,1,
                    1.0,res,1);
        resn = cblas_dnrm2(n, res,1);
        resvec[k]=resn;
        cblas_daxpy(n, alpha, res,1, X,1);
        k++;
    }
    *nbite = k+1;
    free(res);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    int n = *la;
    int ldab = *lab;
    for(int i=0; i<n; i++){
        for(int j=0; j<ldab; j++){
            MB[i * ldab + j] = 0.0;
        }
        MB[i * ldab + 1] = AB[i * ldab + 1];
    }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    int n = *la;
    int ldab = *lab;
    for(int i=0; i<n; i++){
        for(int j=0; j<ldab; j++){
            MB[i * ldab + j] = 0.0;
        }
        MB[i * ldab + 1] = AB[i * ldab + 1];
        MB[i * ldab + 2] = AB[i * ldab + 2];
    }
}

/*void jacobi_1D(double *AB, int *lab, int *la, int *ku, int *kl,double *RHS, double *x, double *tol, int *maxit,double *resvec, int *nbite, double *EX_SOL){
    int n = *la;
    int k = 0;
    double *x_old = (double*) malloc(n*sizeof(double));
    double *res   = (double*) calloc(n,sizeof(double));
    double resnorm=0.0, error=0.0;

    // initial residual
    cblas_dcopy(n, RHS, 1, res, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans,
                n, n, *kl, *ku,
                -1.0, AB, *lab,
                x, 1,
                1.0, res, 1);
    resnorm = cblas_dnrm2(n, res, 1);
    resvec[0] = resnorm;
    if(EX_SOL != NULL){
        double diff=0.0;
        for(int i=0;i<n;i++){
            double d = x[i] - EX_SOL[i];
            diff += d*d;
        }
        error = sqrt(diff);
        printf("Iter=0: residual=%e, error=%e\n", resnorm, error);
    } else {
        printf("Iter=0: residual=%e\n", resnorm);
    }

    while(resnorm > *tol && k < (*maxit)-1) {
        // x_old = x
        cblas_dcopy(n, x, 1, x_old, 1);

        for(int i=0;i<n;i++){
            double Aii = AB[1*n + i]; // main diag row=1
            double sum = RHS[i];

            // lower diag
            if(i>0){
                double A_im1 = AB[0*n + i]; // row=0
                sum -= A_im1 * x_old[i-1];
            }
            // upper diag
            if(i<n-1){
                double A_ip1 = AB[2*n + i]; // row=2
                sum -= A_ip1 * x_old[i+1];
            }
            x[i] = sum / Aii;
        }

        // recompute residual
        cblas_dcopy(n, RHS, 1, res, 1);
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
                    n, n, *kl, *ku,
                    -1.0, AB, *lab,
                    x, 1,
                    1.0, res, 1);
        resnorm = cblas_dnrm2(n, res, 1);
        k++;
        resvec[k] = resnorm;

        if(EX_SOL != NULL){
            double diff=0.0;
            for(int i=0;i<n;i++){
                double dd = x[i] - EX_SOL[i];
                diff += dd*dd;
            }
            error = sqrt(diff);
            printf("Iter=%d: residual=%e, error=%e\n", k, resnorm, error);
        } else {
            printf("Iter=%d: residual=%e\n", k, resnorm);
        }
    }

    *nbite = k;
    free(x_old);
    free(res);
}

void gauss_seidel_1D(double *AB, int *lab, int *la, int *ku, int *kl,double *RHS, double *x, double *tol, int *maxit,double *resvec, int *nbite, double *EX_SOL){
    int n = *la;
    int k = 0;
    double *res = (double*) calloc(n,sizeof(double));
    double resnorm=0.0, error=0.0;

    // initial residual
    cblas_dcopy(n, RHS, 1, res, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans,
                n, n, *kl, *ku,
                -1.0, AB, *lab,
                x, 1,
                1.0, res, 1);
    resnorm = cblas_dnrm2(n, res, 1);
    resvec[0] = resnorm;

    if(EX_SOL != NULL){
        double diff=0.0;
        for(int i=0; i<n; i++){
            double dd = x[i] - EX_SOL[i];
            diff += dd*dd;
        }
        error = sqrt(diff);
        printf("Iter=0: residual=%e, error=%e\n", resnorm, error);
    } else {
        printf("Iter=0: residual=%e\n", resnorm);
    }

    while(resnorm > *tol && k < (*maxit)-1)
    {
        for(int i=0;i<n;i++){
            double Aii = AB[1*n + i]; // main diag
            double sum = RHS[i];

            // lower diag
            if(i>0){
                double A_im1 = AB[0*n + i];
                sum -= A_im1 * x[i-1];
            }
            // upper diag
            if(i<n-1){
                double A_ip1 = AB[2*n + i];
                sum -= A_ip1 * x[i+1];
            }
            x[i] = sum / Aii;
        }

        // new residual
        cblas_dcopy(n, RHS, 1, res, 1);
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
                    n, n, *kl, *ku,
                    -1.0, AB, *lab,
                    x, 1,
                    1.0, res, 1);
        resnorm = cblas_dnrm2(n, res, 1);
        k++;
        resvec[k] = resnorm;

        if(EX_SOL != NULL){
            double diff=0.0;
            for(int i=0; i<n; i++){
                double dd = x[i] - EX_SOL[i];
                diff += dd*dd;
            }
            error = sqrt(diff);
            printf("Iter=%d: residual=%e, error=%e\n", k, resnorm, error);
        } else {
            printf("Iter=%d: residual=%e\n", k, resnorm);
        }
    }

    *nbite = k;
    free(res);
}

void jacobi_solve_Mz(const double* r, double* z, int n)
{
    // for i=0..n-1: z[i] = r[i]/2.0
    for(int i=0; i<n; i++){
        z[i] = 0.5 * r[i];
    }
}

void gs_solve_Mz(const double* MB, int ldab, int n,
                 const double* r, double* z)
{
    // MB[row + col*ldab], row=1 => diag=2, row=2 => -1
    //  (we do not pivot; we trust diag=2 => invertible)
    z[0] = r[0] / MB[1 + 0*ldab]; // = 2
    for(int i=1; i<n; i++){
        double diag = MB[1 + i*ldab]; // row=1 => diag=2
        double low  = MB[2 + i*ldab]; // row=2 => -1
        // z[i] = (r[i] - low*z[i-1]) / diag
        // but for GS => M= (D-E), so Mz = r => z[i] = (r[i] + z[i-1]) / 2 if low=-1
        z[i] = (r[i] + z[i-1]*(-low)) / diag;
        // or (r[i] - (-1)*z[i-1]) /2 => (r[i] + z[i-1]) /2
    }
}
*/

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite,int *NRHS, int *ipiv, int *info){
    int n = *la;
    double *res = (double*) malloc(n*sizeof(double));
    int k=0;
    int ku_m=*ku-1;
    double resn=1.0;

    dgbtrf_(la, la, kl, &ku_m, MB, lab, ipiv, info);

    while((resn>*tol) && (k<(*maxit)-1)){
        // new res
        cblas_dcopy(n, RHS,1,res,1);
        cblas_dgbmv(CblasColMajor, CblasNoTrans,
                    n, n, *kl, *ku,
                    -1.0, AB,*lab,
                    X,1,
                    1.0,res,1);
        resn = cblas_dnrm2(n, res,1);
        resvec[k]=resn;
        dgbtrs_("N", la, kl, &ku_m, NRHS, MB, lab, ipiv, res, la, info, (unsigned long)1);
        cblas_daxpy(n, 1.0, res,1, X,1);
        k++;
    }
    *nbite = k+1;
    free(res);
}

/*void plot_convergence(double *resvec, int nbite) {
    FILE *gnuplot = popen("gnuplot -persist", "w");
    if (gnuplot != NULL) {
        fprintf(gnuplot, "set title 'Convergence History'\n");
        fprintf(gnuplot, "set xlabel 'Iteration'\n");
        fprintf(gnuplot, "set ylabel 'Residual Norm'\n");

        fprintf(gnuplot, "plot '-' with linespoints title 'Residual Norm'\n");

        for (int i = 0; i < nbite; i++) {
            fprintf(gnuplot, "%d %lf\n", i+1, resvec[i]);
        }

        fprintf(gnuplot, "e\n");

        fflush(gnuplot);
        printf("Press Enter to close the plot window...\n");
        getchar();

        pclose(gnuplot);
    } else {
        printf("Error: Could not open gnuplot.\n");
    }
}
*/
void build_poisson1D_CSR(int n, int** csrRowPtr, int** csrColInd, double** csrVal) {
    int nnz = 3*n - 2;
    *csrRowPtr = (int*)   malloc((n+1)*sizeof(int));
    *csrColInd = (int*)   malloc(nnz*sizeof(int));
    *csrVal    = (double*)malloc(nnz*sizeof(double));

    int* rowPtr = *csrRowPtr;
    int* colInd = *csrColInd;
    double* val = *csrVal;

    int idx = 0;
    rowPtr[0] = 0;

    for(int i=0; i<n; i++){
        // if i>0 => A[i,i-1] = -1
        if(i>0) {
            colInd[idx] = i-1;
            val[idx] = -1.0;
            idx++;
        }

        // A[i,i] = 2
        colInd[idx] = i;
        val[idx] = 2.0;
        idx++;

        // If i<n-1 => A[i,i+1] = -1
        if(i < n-1) {
            colInd[idx] = i+1;
            val[idx] = -1.0;
            idx++;
        }

        rowPtr[i+1] = idx;
    }
}

void build_poisson1D_CSC(int n, int** cscColPtr, int** cscRowInd, double** cscVal){
    int nnz = 3*n - 2;
    *cscColPtr = (int*)   malloc((n+1)*sizeof(int));
    *cscRowInd = (int*)   malloc(nnz*sizeof(int));
    *cscVal    = (double*)malloc(nnz*sizeof(double));

    int* colPtr = *cscColPtr;
    int* rowInd = *cscRowInd;
    double* val = *cscVal;

    int idx = 0;
    colPtr[0] = 0;

    for(int j=0; j<n; j++) {
        // if j>0 => A[j-1, j] = -1 => position (j-1)
        if(j>0){
            rowInd[idx] = j-1;
            val[idx] = -1.0;
            idx++;
        }

        // A[j,j] = 2
        rowInd[idx] = j;
        val[idx] = 2.0;
        idx++;

        // If j<n-1 => A[j+1, j] = -1 => position (j+1)
        if(j < n-1) {
            rowInd[idx] = j+1;
            val[idx] = -1.0;
            idx++;
        }

        colPtr[j+1] = idx;
    }
}

void dcsr_mv(int n,const int *csrRowPtr, const int *csrColInd, const double *csrVal,const double *x, double *y) {
    // First clear y to zero
    for(int i=0; i<n; i++){
        y[i] = 0.0;
    }

    for(int i=0; i<n; i++){
        int start = csrRowPtr[i];
        int end   = csrRowPtr[i+1];
        // Iterate over all non-zero elements in this row
        for(int idx = start; idx<end; idx++){
            int j = csrColInd[idx];  // The non-zero element is in the jth column
            double val = csrVal[idx];
            y[i] += val * x[j];
        }
    }
}

void dcsc_mv(int n,const int *cscColPtr, const int *cscRowInd, const double *cscVal,const double *x, double *y) {
    for(int i=0; i<n; i++){
        y[i] = 0.0;
    }

    // For each column j
    for(int j=0; j<n; j++){
        double xj = x[j];
        int start = cscColPtr[j];
        int end   = cscColPtr[j+1];
        for(int idx = start; idx<end; idx++){
            int i = cscRowInd[idx];
            double val = cscVal[idx];
            y[i] += val * xj;
        }
    }
}