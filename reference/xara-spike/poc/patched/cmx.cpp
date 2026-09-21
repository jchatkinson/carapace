// SRC/matrix/routines/cmx.h declares cmx_inv2..cmx_inv6(_v2/_v3) but no
// implementation exists anywhere in this checkout (git history shows the
// header deleted in commit d97b851 "clean up" and reintroduced later
// without its .cpp -- a pre-existing gap in the tree, unrelated to wasm
// porting). This is a generic Gauss-Jordan stand-in, correctness-first
// rather than the closed-form-per-size version the real file presumably
// had, so the POC's Matrix::Invert()/Solve() paths link and run correctly.
#include <cstring>
#include <algorithm>

static int cmx_inv_n(const double *a, double *ainv, int n, int *ok_flag) {
  double m[36 * 2];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) m[i * 2 * n + j] = a[i * n + j];
    for (int j = 0; j < n; j++) m[i * 2 * n + n + j] = (i == j) ? 1.0 : 0.0;
  }
  for (int col = 0; col < n; col++) {
    int piv = col;
    double best = m[col * 2 * n + col] < 0 ? -m[col * 2 * n + col] : m[col * 2 * n + col];
    for (int r = col + 1; r < n; r++) {
      double v = m[r * 2 * n + col] < 0 ? -m[r * 2 * n + col] : m[r * 2 * n + col];
      if (v > best) { best = v; piv = r; }
    }
    if (best == 0.0) { *ok_flag = 0; return -1; }
    if (piv != col)
      for (int k = 0; k < 2 * n; k++) std::swap(m[col * 2 * n + k], m[piv * 2 * n + k]);
    double d = m[col * 2 * n + col];
    for (int k = 0; k < 2 * n; k++) m[col * 2 * n + k] /= d;
    for (int r = 0; r < n; r++) {
      if (r == col) continue;
      double f = m[r * 2 * n + col];
      if (f == 0.0) continue;
      for (int k = 0; k < 2 * n; k++) m[r * 2 * n + k] -= f * m[col * 2 * n + k];
    }
  }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) ainv[i * n + j] = m[i * 2 * n + n + j];
  *ok_flag = 1;
  return 0;
}

extern "C" {
int cmx_inv2(const double *a, double *ainv, int *ok_flag) { return cmx_inv_n(a, ainv, 2, ok_flag); }
int cmx_inv3(const double *a, double *ainv, int *ok_flag) { return cmx_inv_n(a, ainv, 3, ok_flag); }
int cmx_inv4(const double *a, double *ainv, int *ok_flag) { return cmx_inv_n(a, ainv, 4, ok_flag); }
int cmx_inv5(const double *a, double *ainv, int *ok_flag) { return cmx_inv_n(a, ainv, 5, ok_flag); }
int cmx_inv6(const double *a, double *ainv, int *ok_flag) { return cmx_inv_n(a, ainv, 6, ok_flag); }
int cmx_inv6_v2(const double *a, double *ainv, int *ok_flag) { return cmx_inv_n(a, ainv, 6, ok_flag); }
int cmx_inv6_v3(const double *a, double *ainv, int *ok_flag) { return cmx_inv_n(a, ainv, 6, ok_flag); }
}
