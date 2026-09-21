// Stage 1 spike: verify f2c-translated legacy ARPACK (dsaupd/dseupd) against
// a known-eigenvalue problem, compiled as plain C++ (no Fortran toolchain).
//
// Operator: A = diag(1, 2, 3, ..., N)  =>  eigenvalues are exactly 1..N.
// We ask ARPACK for the NEV largest-magnitude eigenvalues and check them
// against the known analytic answer.

#include "f2c.h"
#undef abs
#undef min
#undef max

#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

extern "C" {
int dsaupd_(integer *ido, char *bmat, integer *n, char *which, integer *nev,
            doublereal *tol, doublereal *resid, integer *ncv, doublereal *v,
            integer *ldv, integer *iparam, integer *ipntr, doublereal *workd,
            doublereal *workl, integer *lworkl, integer *info,
            ftnlen bmat_len, ftnlen which_len);

int dseupd_(logical *rvec, char *howmny, logical *select, doublereal *d__,
            doublereal *z__, integer *ldz, doublereal *sigma, char *bmat,
            integer *n, char *which, integer *nev, doublereal *tol,
            doublereal *resid, integer *ncv, doublereal *v, integer *ldv,
            integer *iparam, integer *ipntr, doublereal *workd,
            doublereal *workl, integer *lworkl, integer *info,
            ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len);
}

static void av(int n, const double *x, double *y) {
  for (int i = 0; i < n; i++) y[i] = double(i + 1) * x[i];
}

int main() {
  const integer n = 100;
  const integer nev = 6;
  const integer ncv = 20;
  const integer ldv = n;

  std::vector<doublereal> resid(n), v(ldv * ncv);
  std::vector<doublereal> workd(3 * n);
  integer lworkl = ncv * (ncv + 8);
  std::vector<doublereal> workl(lworkl);
  integer iparam[11] = {0};
  integer ipntr[11] = {0};

  char bmat = 'I';
  char which[3] = "LM";
  doublereal tol = 0.0;
  integer info = 0;
  integer ido = 0;

  iparam[0] = 1;      // exact shifts
  iparam[2] = 300;    // max iterations
  iparam[3] = 1;
  iparam[6] = 1;      // standard eigenproblem, A*x = lambda*x

  int iter = 0;
  for (;;) {
    dsaupd_(&ido, &bmat, (integer *)&n, which, (integer *)&nev, &tol,
            resid.data(), (integer *)&ncv, v.data(), (integer *)&ldv, iparam,
            ipntr, workd.data(), workl.data(), &lworkl, &info, 1L, 2L);
    if (ido == -1 || ido == 1) {
      av(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
      if (++iter > 10000) { printf("FAIL: too many reverse-comm iterations\n"); return 1; }
      continue;
    }
    break;
  }

  if (info < 0) {
    printf("FAIL: dsaupd info=%ld\n", (long)info);
    return 1;
  }

  logical rvec = 0;   // eigenvalues only
  char howmny = 'A';
  std::vector<logical> select(ncv);
  std::vector<doublereal> d(nev);
  doublereal sigma = 0.0;

  dseupd_(&rvec, &howmny, select.data(), d.data(), v.data(), (integer *)&ldv,
          &sigma, &bmat, (integer *)&n, which, (integer *)&nev, &tol,
          resid.data(), (integer *)&ncv, v.data(), (integer *)&ldv, iparam,
          ipntr, workd.data(), workl.data(), &lworkl, &info, 1L, 1L, 2L);

  if (info != 0) {
    printf("FAIL: dseupd info=%ld\n", (long)info);
    return 1;
  }

  // Known answer: largest-magnitude eigenvalues of diag(1..100) are
  // {100, 99, 98, 97, 96, 95} in some order.
  std::vector<double> got(d.begin(), d.end());
  std::sort(got.begin(), got.end());
  double expected[6] = {95, 96, 97, 98, 99, 100};

  printf("ARPACK converged eigenvalues (sorted): ");
  for (auto x : got) printf("%.6f ", x);
  printf("\n");

  bool ok = true;
  for (int i = 0; i < nev; i++) {
    if (std::fabs(got[i] - expected[i]) > 1e-8) ok = false;
  }
  printf(ok ? "PASS\n" : "FAIL: mismatch vs expected {95,96,97,98,99,100}\n");
  return ok ? 0 : 1;
}
