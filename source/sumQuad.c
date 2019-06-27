#include <stdio.h>


int main() {
	long double x, s=0.0;

	do {
	scanf("%LE\n", &x);
	//printf("%40.35LE\n", x);
	s += x;
	} while (!feof(stdin));

	printf("%40.35LE\n", s);
	return 1;
}
