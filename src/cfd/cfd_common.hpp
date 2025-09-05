#ifndef CFD_COMMON_HPP
#define CFD_COMMON_HPP

#include "macros.hpp"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfd {

/**
 * @brief invalid index
 *
 * Used in grid connectivity tables to define blank connection
 */
constexpr size_t INVALID_INDEX = (size_t)-1;

int AAA = 56;
} // namespace cfd

#endif
