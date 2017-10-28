/*
 Copyright (C) 2015 Nilakantha Paudel, Loukas Georgiadis, Giuseppe F. Italiano
 All rights reserved.

 This software may be modified and distributed under the terms
 of the BSD license as followoing.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation  and/or other materials provided with the distribution.

 3. Neither the names of the copyright holders nor the names of any
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
 */

//============================================================================
// Project     : TWO VERTEX STRONGLY CONNECTED COMPONENT
// Name		   : common_system_header.h
// Author      : Nilakantha Paudel [Nilu]
// Version     : 1.0
// Description : Even though, all the header file will have the safeguard.
// However, to keep the code neat and clean, we imported the common system header
// file once in our project and used it as a global.
//============================================================================


#ifndef TWOVSCC_COMMON_COMMON_SYSTEM_HEADER_H_
#define TWOVSCC_COMMON_COMMON_SYSTEM_HEADER_H_

//* Standard System Library on alphabetical order */
#include <algorithm>
#include <assert.h>
#include <cassert>
#include <climits> // for minimum and maximum primitive value
#include <cmath> //for log
#include <cstring> // for memcpy
#include <cstdlib> // for standard library (exit function, division and etc)
#include <fstream> // read write file
#include <iostream> // read file
#include <iomanip> // for set precision (rounding) operation
#include <iterator>
#include <list>
#include <queue> // queue for bfs
#include <set>
#include <map>
#include <string> // for stoi function
#include <sstream>
#include <vector>

// Following for the system timer to sleep the thread and get the random number
#include <ctime>
#include <thread>
#include <chrono>

/*-----------------------------------------------------------------------*
 | invoke common std parameter rather then using the whole std namespace |
 | -- an alphabetical order												 |
 *-----------------------------------------------------------------------*/

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::ofstream;
using std::exception;
using std::min;

// to fixed the decimal precision
//using std::fixed;
// map function
//using std::map;
// to read and write the file
//using std::ofstream;
// stream template
//using std::ostream;
//pointer to the functional object
//using std::ptr_fun;
// to create the pair
//using std::pair;
// to create the queue
//using std::queue;
// to set up the decimal percison
//using std::setprecision;
// to swap the values between two variables
//using std::swap;



#endif /* TWOVSCC_COMMON_COMMON_SYSTEM_HEADER_H_ */
