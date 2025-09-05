#ifndef CFD24_MACROS_HPP
#define CFD24_MACROS_HPP

#include <stdexcept>
#include <sstream>

#if !defined(__PRETTY_FUNCTION__) && !defined(__GNUC__)
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#define _THROW_NOT_IMP_ \
{\
	std::ostringstream oss; \
	oss << "NOT IMPLEMENTED method called." << std::endl; \
	oss << __PRETTY_FUNCTION__ << std::endl; \
	oss << "At" << std::endl; \
	oss << __FILE__ << ": " << __LINE__ << std::endl; \
	throw std::runtime_error(oss.str()); \
}

#define _THROW_UNREACHABLE_ \
{\
	std::ostringstream oss; \
	oss << "Unreachable code." << std::endl; \
	oss << __PRETTY_FUNCTION__ << std::endl; \
	oss << "At" << std::endl; \
	oss << __FILE__ << ": " << __LINE__ << std::endl; \
	throw std::runtime_error(oss.str()); \
}

#define _THROW_INTERNAL_ERROR_ \
{\
	std::ostringstream oss; \
	oss << "Internal error." << std::endl; \
	oss << __PRETTY_FUNCTION__ << std::endl; \
	oss << "At" << std::endl; \
	oss << __FILE__ << ": " << __LINE__ << std::endl; \
	throw std::runtime_error(oss.str()); \
}

#endif
