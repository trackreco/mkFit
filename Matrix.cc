#include "Matrix.h"


std::default_random_engine            g_gen(0xbeef0133);
std::normal_distribution<float>       g_gaus(0.0, 1.0);
std::uniform_real_distribution<float> g_unif(0.0, 1.0);
