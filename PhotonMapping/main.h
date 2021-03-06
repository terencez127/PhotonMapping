//
//  main.h
//  PhotonMapping
//
//  Copyright (c) 2014 Terence. All rights reserved.
//

#ifndef PhotonMapping_main_h
#define PhotonMapping_main_h

#include <cmath>
#include <cassert>

namespace pm {
    /// Color class
    struct color {
        unsigned char rgba[4];
        color () {}
        color (double r, double g, double b, double a = 255);
        color (double gray, double alpha = 255);
    };

    void rect (int x, int y, int a, int b);

}

void reset();
void display_scene();
void display_photon();

#endif
