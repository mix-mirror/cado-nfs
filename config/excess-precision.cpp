/* Does the compiler evaluate a double as a double?
 *
 * On i386 with x87 arithmetic -- which is what gcc selects by default
 * for -m32, since -mfpmath defaults to 387 there -- intermediate results
 * are kept in 80-bit registers and are not rounded back to double at
 * every step. Two expressions that are mathematically equal can then
 * compare unequal, and code that decides anything on such a comparison
 * does not give the same answer as it does elsewhere.
 *
 * C says which of the two we are in: FLT_EVAL_METHOD is 0 when every
 * operation is evaluated at its own type, and 1 or 2 when it is
 * evaluated at a wider one.
 */

#include <cfloat>

#if !defined(FLT_EVAL_METHOD) || FLT_EVAL_METHOD != 0
#error "double arithmetic carries excess precision"
#endif

int main()
{
    return 0;
}
