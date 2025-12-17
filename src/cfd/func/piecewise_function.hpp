#ifndef __CFD_FUNC_PIECEWISE_FUNCTION_HPP__
#define __CFD_FUNC_PIECEWISE_FUNCTION_HPP__

#include "cfd/func/i_function.hpp"

namespace cfd {

template<typename F>
class PiecewiseFunction : public F {
public:
};

} // namespace cfd

#endif
