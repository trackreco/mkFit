//#define DEBUG
#ifdef DEBUG
#define dprint(x) if (debug) std::cout << x << std::endl
#define dcall(x)  if (debug) { x; }
#else
#define dprint(x) (void(0))
#define dcall(x) (void(0))
#endif
