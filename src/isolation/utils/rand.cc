#include "../utils/utils.hpp"
#include <stdint.h>
#include <ctime>

Rand randgen(walltime() * 1e9);

// Ranq1, from "Numerical Recipes 3rd edition,"
// Press, Tenkolsky, Vetterling and Flannery, 2007.

static const uint64_t Mul = 2685821657736338717LL;
static const uint64_t Vini = 4101842887655102017LL;
static const double Fl = 5.42101086242752217e-20;

Rand::Rand(uint64_t seed) : theseed(seed) {
	v = Vini ^ seed;
	v = bits();
}

uint64_t Rand::bits() {
	v ^= v >> 21;
	v ^= v << 35;
	v ^= v >> 4;
	v *= Mul;
	return v;
}

int64_t Rand::integer(long min, long max) {
	if (min == max)
		return min;
	
	return bits() % ((max+1) - min) + min;
}

double Rand::real() {
	return bits() * Fl;
}