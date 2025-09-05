#ifndef CFD_COMMON_HPP
#define CFD_COMMON_HPP

#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <memory>
#include <iostream>
#include <sstream>
#include "macros.hpp"

namespace cfd{

/**
 * @brief invalid index
 *
 * Used in grid connectivity tables to define blank connection
 */
constexpr size_t INVALID_INDEX = (size_t)-1;

}

#endif
