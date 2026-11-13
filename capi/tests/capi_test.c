/* Minimal C test of the Stratonix C API: validate, build, evaluate, error handling. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "stratonix_c.h"

static int failures = 0;
#define EXPECT(cond, msg)                                     \
  do {                                                        \
    if (!(cond)) {                                            \
      printf("FAIL %s (%s) last error: %s\n", msg, #cond, stx_last_error()); \
      ++failures;                                             \
    } else {                                                  \
      printf("PASS %s\n", msg);                               \
    }                                                         \
  } while (0)

int main(int argc, char** argv) {
  if (argc < 3) {
    printf("usage: capi_test <valid description> <invalid description>\n");
    return 2;
  }
  EXPECT(stx_capi_version() >= 100, "C API version");
  EXPECT(strlen(stx_library_version()) > 0, "library version string");
  EXPECT(stx_validate(argv[1]) == STX_OK, "valid description validates");
  EXPECT(stx_validate(argv[2]) == STX_E_VALIDATION, "invalid description is refused");
  EXPECT(strlen(stx_last_error()) > 0, "error message is set");
  EXPECT(stx_validate(NULL) == STX_E_INVALID_ARGUMENT, "NULL path is refused");

  stx_model* m = NULL;
  EXPECT(stx_model_build(argv[1], &m) == STX_OK && m != NULL, "model builds");
  if (m) {
    int nh = 0, nu = 0, nf = 0;
    EXPECT(stx_model_counts(m, &nh, &nu, &nf) == STX_OK && nh > 0 && nu == nh + 1, "counts");
    EXPECT(stx_model_horizon_name(m, 0) != NULL, "horizon name");
    EXPECT(stx_model_horizon_name(m, 999) == NULL, "out-of-range horizon name is NULL");
    double xyz[6] = {1000.0, 1000.0, -10.0, 1000.0, 1000.0, -2400.0};
    int units[2] = {-1, -1};
    double* scal = (double*)0;
    double s[64];
    scal = s;
    EXPECT(stx_model_evaluate(m, 2, xyz, units, NULL, scal) == STX_OK, "evaluate");
    EXPECT(units[0] == 0 && units[1] == nh, "shallow point above top, deep point below base");
    EXPECT(s[0] > 0 && s[1 * 1 + 0] < 0, "scalar signs (horizon-major)");
    stx_model_free(m);
  }
  printf("%s\n", failures ? "C API TEST FAILED" : "C API TEST PASSED");
  return failures ? 1 : 0;
}
