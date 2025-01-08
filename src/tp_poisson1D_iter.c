/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <time.h>

#define ALPHA 0
#define JAC 1
#define GS 2

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
    int nbpoints, la;
    int ku, kl, lab, kv;
    int *ipiv=NULL;
    int info;
    int NRHS;
    int IMPLEM = 0;
    double T0, T1;
    double *RHS, *SOL, *EX_SOL, *X;
    double *AB;
    double *MB;

    double relres;

    double opt_alpha;

    if (argc == 2) {
        IMPLEM = atoi(argv[1]);
    } else if (argc > 2) {
        perror("Application takes at most one argument");
        exit(1);
    }

    /* Size of the problem */
    NRHS=1;
    nbpoints=12;
    la=nbpoints-2;

    /* Dirichlet Boundary conditions */
    T0=5.0;
    T1=20.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS=(double *) malloc(sizeof(double)*la);
    SOL=(double *) calloc(la, sizeof(double));
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

    /* Setup the Poisson 1D problem */
    /* General Band Storage */
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    write_vec(RHS, &la, "RHS.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(X, &la, "X_grid.dat");

    kv = 0;
    ku = 1;
    kl = 1;
    lab = kv + kl + ku + 1;
    AB = (double *) malloc(sizeof(double)*lab*la);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    /* uncomment the following to check matrix A */
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

    /********************************************/
    /* Solution (Richardson with optimal alpha) */

    /* Computation of optimum alpha */
    opt_alpha = richardson_alpha_opt(&la);
    printf("Optimal alpha for simple Richardson iteration is : %lf",opt_alpha);

    /* Solve */
    double tol=1e-3;
    int maxit=1000;
    double *resvec;
    int nbite=0;

    resvec=(double *) calloc(maxit, sizeof(double));

    // Performance Testing: Timing
    clock_t start, end;
    double time_spent = 0.0;

    /* Solve with Richardson alpha */
    if (IMPLEM == ALPHA) {
        start = clock();
        richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
        end = clock();
        time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("[Richardson alpha] CPU time = %f s\n", time_spent);
        printf("[Richardson alpha] Iterations = %d\n", nbite);
        write_vec(resvec, &nbite, "ALPHA_RESVEC.dat");
    }
    /*else if (IMPLEM==JAC){
        start = clock();
        jacobi_1D(AB, &lab, &la, &ku, &kl,RHS, SOL, &tol, &maxit,resvec, &nbite, EX_SOL);
        end = clock();
        time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("[Jacobi] CPU time = %f s\n", time_spent);
        printf("[Jacobi] Iterations = %d\n", nbite);
    }
    else if (IMPLEM == GS) {
        start = clock();
        gauss_seidel_1D(AB, &lab, &la, &ku, &kl,RHS, SOL, &tol, &maxit,resvec, &nbite, EX_SOL);
        end = clock();
        time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("[Gauss-Seidel] CPU time = %f s\n", time_spent);
        printf("[Gauss-Seidel] Iterations = %d\n", nbite);
    }*/

    kv = 1;
    ku = 1;
    kl = 1;
    MB = (double *) malloc(sizeof(double) * lab * la);
    ipiv = (int *) malloc(la * sizeof(int));
    /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
    if (IMPLEM == JAC) {
        extract_MB_jacobi_tridiag(AB, MB, &lab,&la,&ku,&kl,&kv);
        write_GB_operator_colMajor_poisson1D(MB,&lab,&la,"MB_jac.dat");
        start=clock();
        print_GB_operator_colMajor_poisson1D(MB, lab, la, "MB BEFORE calling richardson_MB");
        richardson_MB(AB,RHS,SOL,MB,&lab,&la,&ku,&kl,&tol,&maxit,resvec,&nbite,&NRHS,ipiv,&info);
        print_GB_operator_colMajor_poisson1D(MB, lab, la, "MB BEFORE calling richardson_MB");
        end=clock();
        time_spent=(double)(end-start)/CLOCKS_PER_SEC;
        printf("[Jacobi => M=D => richardson_MB] CPU time= %f s\n", time_spent);
        printf("Iterations= %d\n", nbite);
        write_vec(resvec, &nbite, "JAC_RESVEC.dat");

    } else if (IMPLEM == GS) {
        extract_MB_gauss_seidel_tridiag(AB, MB, &lab,&la,&ku,&kl,&kv);
        write_GB_operator_colMajor_poisson1D(MB,&lab,&la,"MB_gs.dat");
        start=clock();
        print_GB_operator_colMajor_poisson1D(MB, lab, la, "MB BEFORE calling richardson_MB");
        richardson_MB(AB,RHS,SOL,MB,&lab,&la,&ku,&kl,&tol,&maxit,resvec,&nbite,&NRHS,ipiv,&info);
        print_GB_operator_colMajor_poisson1D(MB, lab, la, "MB BEFORE calling richardson_MB");
        end=clock();
        time_spent=(double)(end-start)/CLOCKS_PER_SEC;
        printf("[GaussSeidel => M=(D-E) => richardson_MB] CPU time= %f s\n", time_spent);
        printf("Iterations= %d\n", nbite);
        write_vec(resvec, &nbite, "GS_RESVEC.dat");
    }

    /* Write solution */
    write_vec(SOL, &la, "SOL.dat");

    /* Relative forward error */
    relres = relative_forward_error(SOL, EX_SOL, &la);
    printf("\nThe relative forward error is relres = %e\n",relres);

    //plot_convergence(resvec,nbite);

    /*free(resvec);
    free(RHS);
    free(SOL);
    free(EX_SOL);
    free(X);
    free(AB);
    free(MB);
    */
    printf("\n\n--------- End -----------\n");
}