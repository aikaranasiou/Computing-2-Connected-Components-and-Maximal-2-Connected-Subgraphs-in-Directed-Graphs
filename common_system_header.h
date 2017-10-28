


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
