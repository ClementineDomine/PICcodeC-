//
//  binary.hpp
//  Danielpenningtrap
//
//  Created by Clementine on 07.03.20.
//  Copyright © 2020 Clementine. All rights reserved.
//

#ifndef binary_hpp
#define binary_hpp

#include <stdio.h>

#include <fstream.h>

char buffer[100];
ofstream myFile ("data.bin", ios::out | ios::binary);
myFile.write (buffer, 100);

#endif /* binary_hpp */
