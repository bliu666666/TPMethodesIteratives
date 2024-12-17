/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include <time.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X,*Y;
  double **AAB;
  double *AB;

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  Y=(double *) malloc(sizeof(double)*la); // Utilisé pour stocker le résultat du produit matrice-vecteur

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    // 打印 AB 看是否主对角线全为 2
    printf("---- Print AB (colMajor) after set_GB_operator_colMajor_poisson1D ----\n");
    for (int row=0; row<lab; row++){
        for (int col=0; col<la; col++){
            printf("%+6.2f ", AB[row*la + col]);
        }
        printf("\n");
    }

    printf("\nMain diagonal:\n");
    for(int i=0; i<la; i++){
        printf(" AB[kv*la + i] = AB[%d] = %+6.2f\n", i, AB[kv*la + i]);
    }
    printf("-----------------------------------------------------\n");
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  // Multiplication matrice-vecteur à l'aide de la fonction dgbmv
  printf("Performing matrix-vector multiplication with dgbmv...\n");
  cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1.0,AB,lab,EX_SOL,1,0.0,Y,1);

  write_vec(Y, &la, "MatrixVectorResult.dat"); // Écrire le produit matrice-vecteur dans un fichier

  // Résultats de la vérification (erreur relative)
  double error = relative_forward_error(Y, RHS, &la);
  printf("Relative forward error for dgbmv = %e\n", error);

  if (error < 1e-10) {
      printf("Matrix-vector multiplication verified successfully.\n");
  } else {
      printf("Matrix-vector multiplication verification failed. Error = %e\n", error);
  }

  printf("Solution with LAPACK\n");
  ipiv = (int *)calloc(la,sizeof(int));
  clock_t start, end;
  double lu_time, solve_time;
  /* LU Factorization */
  if (IMPLEM == TRF) {
      start = clock();
      dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
      end = clock();
      lu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
      if (info==0)
          printf("LU Factorization (Custom) completed successfully in %f seconds.\n", lu_time);
      else
          printf("Error during custom LU Factorization. INFO = %d\n", info);
  }

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if (IMPLEM == TRI) {
      start = clock();
      info = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
      end = clock();
      lu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
      if (info==0)
          printf("LU Factorization (Custom) completed successfully in %f seconds.\n", lu_time);
      else
          printf("Error during custom LU Factorization. INFO = %d\n", info);
  }

  if (IMPLEM == TRI || IMPLEM == TRF){
    //Solution (Triangular)
    if (info==0){
        start = clock();
        LAPACK_dgbtrs("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        end = clock();
        solve_time = ((double)(end - start)) / CLOCKS_PER_SEC;
        if (info == 0) {
            printf("Solution computed successfully in %f seconds.\n", solve_time);
        } else {
            printf("Error during triangular solve. INFO = %d\n", info);
        }
    }else{
      printf("\n INFO = %d\n",info);
    }
  }


  /* It can also be solved with dgbsv */
  if (IMPLEM == SV) {
      // TODO : use dgbsv
      start = clock();
      dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      end = clock();
      lu_time = ((double)(end - start)) / CLOCKS_PER_SEC;

      if (info == 0) {
          printf("Direct solution using dgbsv completed successfully in %f seconds.\n", lu_time);
      } else {
          printf("Error during direct solve. INFO = %d\n", info);
      }
  }

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  /*if (AB) free(AB);
  if (X) free(X);
  if (EX_SOL) free(EX_SOL);
  if (RHS) free(RHS);
  //if (Y) free(Y);
  //if (ipiv) free(ipiv);
   */
  printf("\n\n--------- End -----------\n");
}
