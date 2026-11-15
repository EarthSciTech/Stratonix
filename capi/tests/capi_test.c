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
  if (argc < 5) {
    printf("usage: capi_test <valid description> <invalid description> <faulted description> <new temp folder>\n");
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
  EXPECT(strcmp(stx_status_name(STX_E_CHECK_FAILED), "STX_E_CHECK_FAILED") == 0, "status names");

  /* Faulted model: names, surfaces with gaps, per-block surface, export, package. */
  m = NULL;
  EXPECT(stx_model_build(argv[3], &m) == STX_OK && m != NULL, "faulted model builds");
  if (m) {
    int nh = 0, nu = 0, nf = 0;
    stx_model_counts(m, &nh, &nu, &nf);
    EXPECT(nf == 2 && strcmp(stx_model_fault_name(m, 1), "F2") == 0, "fault names");
    EXPECT(stx_model_fault_name(m, 2) == NULL, "out-of-range fault name is NULL");
    static double z[161 * 121];
    EXPECT(stx_model_horizon_grid(m, 1, 0.0, 0.0, 50.0, 50.0, 161, 121, -1, z) == STX_OK, "faulted horizon grid");
    int gaps = 0, defined = 0;
    for (int k = 0; k < 161 * 121; ++k) {
      if (isnan(z[k])) ++gaps;
      else ++defined;
    }
    EXPECT(gaps > 0 && defined > 19000, "fault gaps present, the rest defined");
    EXPECT(stx_model_horizon_grid(m, 1, 0.0, 0.0, 50.0, 50.0, 161, 121, 1, z) == STX_OK && !isnan(z[0]), "per-block grid has no gaps");
    EXPECT(stx_model_horizon_grid(m, 99, 0.0, 0.0, 50.0, 50.0, 2, 2, -1, z) == STX_E_INVALID_ARGUMENT, "bad horizon index refused");
    stx_model_free(m);
  }
  char out[1024];
  snprintf(out, sizeof out, "%s/capi_export", argv[4]);
  EXPECT(stx_export(argv[3], out, "irap_ascii,well_tops") == STX_OK, "export");
  EXPECT(stx_export(argv[3], out, "irap_ascii") == STX_E_IO, "export into a non-empty folder is refused");
  snprintf(out, sizeof out, "%s/capi_package", argv[4]);
  EXPECT(stx_package(argv[3], out, "{\"formats\": [\"csv\"], \"figures\": false, \"zip\": true}") == STX_OK, "package");
  snprintf(out, sizeof out, "%s/capi_package_bad", argv[4]);
  EXPECT(stx_package(argv[3], out, "{\"banned_patterns\": [\"(unclosed\"]}") == STX_E_VALIDATION, "invalid package option refused");
  /* C API 1.2: a parameter set before building (STX-TN-13 section 7) */
  EXPECT(stx_capi_version() >= 102, "C API 1.2");
  stx_description* d = NULL;
  EXPECT(stx_description_load(argv[3], &d) == STX_OK && d != NULL, "description loads");
  if (d) {
    double v = 0;
    char unit[16] = {0};
    EXPECT(stx_description_get_parameter(d, "faults.F1.throw.max", &v, unit, sizeof unit) == STX_OK && fabs(v - 40.0) < 1e-9 && strcmp(unit, "m") == 0,
           "read faults.F1.throw.max (40 m)");
    EXPECT(stx_description_set_parameter(d, "faults.F1.throw.max", 0.06, "km") == STX_OK, "set faults.F1.throw.max to 0.06 km");
    EXPECT(stx_description_get_parameter(d, "faults.F1.throw.max", &v, NULL, 0) == STX_OK && fabs(v - 60.0) < 1e-9, "read back 60 m");
    EXPECT(stx_description_set_parameter(d, "faults.F1.throw.max", 5.0, "deg") == STX_E_INVALID_ARGUMENT, "a unit of another dimension is refused");
    EXPECT(stx_description_set_parameter(d, "faults.F9.throw.max", 5.0, NULL) == STX_E_INVALID_ARGUMENT, "an unknown path is refused");
    stx_model* m2 = NULL;
    EXPECT(stx_model_build_from(d, &m2) == STX_OK && m2 != NULL, "model builds from the edited description");
    stx_model_free(m2);
    stx_description_free(d);
  }
  printf("%s\n", failures ? "C API TEST FAILED" : "C API TEST PASSED");
  return failures ? 1 : 0;
}
