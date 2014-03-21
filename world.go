package main

// #cgo pkg-config: world
// #include <world/star.h>
import "C"

import (
	"fmt"
)

func main() {
	fmt.Println("Hello world:)")
	a := C.GetFFTSizeForStar(C.int(44100))
	fmt.Println(int(a))
}
