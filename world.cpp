#include <iostream>
#include <world/star.h>

int main() { 
  int a = GetFFTSizeForStar(44100);
  std::cout << "Hello world" << std::endl
	    << a << std::endl; 
  return 0;
}
