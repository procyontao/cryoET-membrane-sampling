#include <cstdio>
#include <cstdarg>

void wprint(const char fmt[], ...) {
	va_list args;
	va_start (args, fmt);
	/*     fmt = va_arg (args, char *); */

	vprintf (fmt, args);

	va_end (args);
}
