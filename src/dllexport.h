//-----------------------------------------------------------------------------
// Ryuichi YAMAMOTO (zryuichi@gmail.com)
//-----------------------------------------------------------------------------

#ifndef WORLD_DLLEXPORT_H_
#define WORLD_DLLEXPORT_H_

#ifdef _WIN32
#    define DLLEXPORT __declspec(dllexport)
#  else
#    define DLLEXPORT
#endif

#endif  // WORLD_DLLEXPORT_H_
