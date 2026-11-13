/* Small example of the Stratonix C API: build a model, read its horizons and faults, evaluate a
 * vertical profile and print the unit at each depth. Build: link against stratonix_c.
 * Usage: stratonix_example <description.yaml> */
#include <stdio.h>
#include <stdlib.h>

#include "stratonix_c.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s <description.yaml>\n", argv[0]);
    return 2;
  }
  printf("Stratonix %s, C API %d\n", stx_library_version(), stx_capi_version());
  stx_model* m = NULL;
  stx_status s = stx_model_build(argv[1], &m);
  if (s != STX_OK) {
    fprintf(stderr, "build failed: %s: %s\n", stx_status_name(s), stx_last_error());
    return 1;
  }
  int nh, nu, nf;
  stx_model_counts(m, &nh, &nu, &nf);
  printf("%d horizons, %d units, %d faults\n", nh, nu, nf);
  for (int f = 0; f < nf; ++f) printf("  fault %d: %s\n", f + 1, stx_model_fault_name(m, f));
  /* A vertical profile at (4000, 3000): 11 points from 600 to 1600 m below the datum. */
  enum { N = 11 };
  double xyz[3 * N];
  int32_t unit[N];
  for (int i = 0; i < N; ++i) {
    xyz[3 * i] = 4000.0;
    xyz[3 * i + 1] = 3000.0;
    xyz[3 * i + 2] = -600.0 - 100.0 * i;
  }
  s = stx_model_evaluate(m, N, xyz, unit, NULL, NULL);
  if (s != STX_OK) {
    fprintf(stderr, "evaluate failed: %s\n", stx_last_error());
    stx_model_free(m);
    return 1;
  }
  for (int i = 0; i < N; ++i) {
    const char* below = unit[i] > 0 ? stx_model_horizon_name(m, unit[i] - 1) : "(above the top horizon)";
    printf("  depth %6.0f m: unit %d, below %s\n", -xyz[3 * i + 2], unit[i], below);
  }
  stx_model_free(m);
  return 0;
}
