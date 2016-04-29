#ifndef _debug_
#define _debug_
#ifdef DEBUG
#include <mutex>

#define dmutex_guard std::lock_guard<std::mutex> dlock(debug_mutex)
#define dprint(x) if (debug) { dmutex_guard; std::cout << x << std::endl; }
#define dcall(x)  if (debug) { dmutex_guard; x; }
#define dprintf(...) if (debug) { dmutex_guard; printf(__VA_ARGS__); }

namespace { 
  bool debug = true; // default, can be overridden locally
  std::mutex debug_mutex;
}

static void print(const TrackState& s)
{
  std::cout << "x:  "  << s.parameters[0] 
            << " y:  " << s.parameters[1]
            << " z:  " << s.parameters[2] << std::endl
            << "px: "  << s.parameters[3]
            << " py: " << s.parameters[4]
            << " pz: " << s.parameters[5] << std::endl
            << "valid: " << s.valid << " errors: " << std::endl;
  dumpMatrix(s.errors);
  std::cout << std::endl;
}

static void print(std::string label, int itrack, const Track& trk)
{
  std::cout << std::endl << label << ": " << itrack << " hits: " << trk.nFoundHits() << " State" << std::endl;
  print(trk.state());
}

static void print(std::string label, const TrackState& s)
{
  std::cout << label << std::endl;
  print(s);
}

static void print(std::string label, const MeasurementState& s)
{
  std::cout << label << std::endl;
  std::cout << "x: "  << s.parameters()[0] 
            << " y: " << s.parameters()[1]
            << " z: " << s.parameters()[2] << std::endl
            << "errors: " << std::endl;
  dumpMatrix(s.errors());
  std::cout << std::endl;
}

#else
#define dprint(x) (void(0))
#define dcall(x) (void(0))
#define dprintf(...) (void(0))
#endif
#endif
