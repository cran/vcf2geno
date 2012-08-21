#include "R.h"

void bz_internal_error ( int errcode ) {
  REprintf("bz_internal_error: code = %d\n", errcode);
}
