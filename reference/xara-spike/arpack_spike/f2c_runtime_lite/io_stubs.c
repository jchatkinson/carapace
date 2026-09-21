/* Minimal stand-ins for f2c's formatted-write runtime (s_wsfe/do_fio/e_wsfe).
 * ARPACK only reaches these through its msglvl-gated debug print calls
 * (ivout/dvout), which are dead code at the default verbosity used in
 * production. Real support would mean porting f2c's FORMAT interpreter
 * (fio.h/fmt.h and friends) -- out of scope for a spike whose only job is
 * to prove the numerics translate correctly.
 */
#include "f2c.h"

/* f2c's own generated call sites declare s_copy as returning int (the
 * historical "Subroutine" convention) even though libf2c defines it void;
 * wasm-ld enforces exact signature matches (traps on call otherwise), so
 * this wrapper matches what the f2c-generated callers actually expect. */
int s_copy(register char *a, register char *b, ftnlen la, ftnlen lb) {
	register char *aend, *bend;
	aend = a + la;
	if (la <= lb) {
		if (a <= b || a >= b + la)
			while (a < aend) *a++ = *b++;
		else
			for (b += la; a < aend; ) *--aend = *--b;
	} else {
		bend = b + lb;
		if (a <= b || a >= bend)
			while (b < bend) *a++ = *b++;
		else {
			a += lb;
			while (b < bend) *--a = *--bend;
			a += lb;
		}
		while (a < aend) *a++ = ' ';
	}
	return 0;
}

int s_wsfe(cilist *a) { (void)a; return 0; }
int e_wsfe(void) { return 0; }
int do_fio(ftnint *number, char *ptr, ftnlen len) { (void)number; (void)ptr; (void)len; return 0; }
