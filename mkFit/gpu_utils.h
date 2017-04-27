#ifndef GPU_UTILS_H
#define GPU_UTILS_H 

#include <cub/util_debug.cuh>
#include <nvToolsExt.h>

#include <cstdint>
#include <unordered_map>

#define cudaCheckError()               \
  do {                                 \
    cudaError_t e=cudaGetLastError();  \
    CubDebugExit(e);                  \
  } while(0)  

#define cudaCheckErrorSync()           \
  do {                                 \
    cudaDeviceSynchronize();           \
    cudaCheckError();                  \
  } while(0)

// CUDA specific:
// Maximum number of blocks in the X direction of the thread grid.
constexpr int max_blocks_x = INT32_MAX;

// The first call to a CUDA API function takes the initialization hit.
void separate_first_call_for_meaningful_profiling_numbers();

void sync_gpu();


class Tracer {
 public:
  Tracer(const char* name, unsigned int color) {
    nvtxEventAttributes_t eventAttrib = {0};
    eventAttrib.version = NVTX_VERSION;
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.colorType = NVTX_COLOR_ARGB;
    eventAttrib.color = color;
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii = name;
    nvtxRangePushEx(&eventAttrib);
  }
  ~Tracer() {
    nvtxRangePop();
  }
};

static std::unordered_map<std::string, unsigned int> tracer_colors = {
    {"Blue",            0x0000FF},
    {"Orange",          0xFFA500},
    {"Red",             0xFF0000},
    {"Pink",            0xFFC0CB},
    {"Green",           0x008000},
    {"PaleGreen",       0x98FB98},
    {"Chocolate",       0xD2691E}
};


#endif /* ifndef GPU_UTILS_H */
