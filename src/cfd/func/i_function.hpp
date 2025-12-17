#ifndef __CFD_FUNC_I_FUNCTION_HPP__
#define __CFD_FUNC_I_FUNCTION_HPP__

namespace cfd {

class IFunction {
public:
    virtual ~IFunction() = default;
    virtual double operator()(double x) const;
};

class IFunctionC1 : public IFunction {
public:
    virtual ~IFunctionC1() = default;
    virtual double dx(double x) const;
};

class IFunctionC2 : public IFunctionC1 {
public:
    virtual ~IFunctionC2() = default;
    virtual double d2x(double x) const;
};

} // namespace cfd

#endif
