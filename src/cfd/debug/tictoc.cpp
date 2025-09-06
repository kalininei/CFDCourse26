#include "tictoc.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace cfd::dbg;

namespace {

/**
 * @brief Execution timer
 *
 * Represents framework for user time measurement.
 *
 * This class can be used as follows: <br>
 *
 * @code{cpp}
 * TicToc timer(); // declare and start timer with default name
 * // ....
 * timer.toc();    // stop timer
 * // ...
 * timer.tic();    // continue timer;
 * // ...
 * timer.fintoc(); // stop and report resulting duration
 * @endcode
 *
 * Also static interface is presented. The below code will result in
 * printing report for main procedure and two subprocedures
 * @code{cpp}
 *
 * Tic("MainProcedure"); // start timer with id = MainProcedure
 * Tic("SubProcedure1"); // start another timer with id = SubProcedure1
 *
 * // .....
 *
 * Toc("SubProcedure1"); // stop SubProcedure1 timer
 * Tic("SubProcedure1"); // start timer with id = SubProcedure2
 *
 * // .....
 *
 * Toc("SubProcedure2"); // stop SubProcedure2 timer
 * Toc("MainProcedure"); // stop MainProcedure timer
 *
 * FinReport();
 * @endcode
 *
 * All static timers exist in global context so they can be called
 * and stopped from different places.
 */
class TicToc {
public:
    /// create time with defined name.
    TicToc(std::string name = "Time duration", bool start = true);

    /// dtor
    ~TicToc() = default;

    /// resets timer. Doesn't restarts it.
    void init();
    /// start timer
    void tic();
    /// pause timer
    void toc();
    /// print report to std::cout
    void report() const;

    /// stop and report to std::cout
    void fintoc() {
        toc();
        report();
    }

    /// get elapsed time in seconds
    double elapsed() const;

private:
    using hr_clock_t = std::chrono::high_resolution_clock;
    using time_point_t = hr_clock_t::time_point;
    using duration_t = std::chrono::duration<double>;

    std::string name_;
    bool is_working_;
    time_point_t tp_;
    duration_t dur_;
};

struct Timers {
    std::map<std::string, TicToc> data;

    ~Timers() {
        FinReport();
    }

    bool has(std::string s) {
        return data.find(s) != data.end();
    }

    TicToc& get(std::string s) {
        auto fnd = data.find(s);
        if (fnd == data.end()) {
            data[s] = TicToc(s, 0);
        }
        return data[s];
    }

    std::vector<std::string> keys() {
        std::vector<std::string> ret;
        for (auto& v : data) {
            ret.push_back(v.first);
        }
        return ret;
    }

    std::vector<std::string> keys_sorted_by_elapsed() {
        std::vector<std::string> ret = keys();
        std::sort(ret.begin(), ret.end(), [this](const std::string& key1, const std::string& key2) -> bool {
            return get(key1).elapsed() > get(key2).elapsed();
        });
        return ret;
    }

    void erase(std::string s) {
        auto fnd = data.find(s);
        if (fnd != data.end())
            data.erase(fnd);
    }
};

Timers alltimers_;

} // namespace

void cfd::dbg::Tic(std::string s) {
    if (s.size() == 0) {
        for (int i = 0; i < 99999; ++i) {
            std::string nm = "Timer" + std::to_string(i);
            if (!alltimers_.has(nm))
                return Tic(nm);
        }
    } else {
        auto& tm = alltimers_.get(s);
        tm.tic();
    }
}

void cfd::dbg::Tic1(std::string s) {
    cfd::dbg::Toc();
    cfd::dbg::Tic(s);
}

void cfd::dbg::Toc(std::string s) {
    if (s.size() == 0) {
        for (auto k : alltimers_.keys())
            Toc(k);
    } else {
        if (alltimers_.has(s)) {
            alltimers_.get(s).toc();
        }
    }
}

void cfd::dbg::Report(std::string s) {
    if (s.size() == 0) {
        for (auto k : alltimers_.keys())
            Report(k);
    } else {
        if (alltimers_.has(s)) {
            alltimers_.get(s).report();
        }
    }
}

void cfd::dbg::FinReport(std::string s) {
    if (s.size() == 0) {
        Toc();
        for (auto k : alltimers_.keys_sorted_by_elapsed()) {
            FinReport(k);
        }
    } else {
        if (alltimers_.has(s)) {
            alltimers_.get(s).fintoc();
        }
        alltimers_.erase(s);
    }
}

TicToc::TicToc(std::string name, bool start) : name_(name), is_working_(start), dur_(duration_t::zero()) {
    if (start)
        tic();
}

void TicToc::init() {
    dur_ = duration_t::zero();
    tp_ = hr_clock_t::now();
}

void TicToc::tic() {
    if (!is_working_) {
        is_working_ = true;
        tp_ = hr_clock_t::now();
    }
}

void TicToc::toc() {
    if (is_working_) {
        is_working_ = false;
        dur_ += std::chrono::duration_cast<duration_t>(hr_clock_t::now() - tp_);
    }
}

void TicToc::report() const {
    std::ostringstream oss;
    oss << std::setw(10) << name_ << ":  " << std::setprecision(3) << std::setw(5) << std::left << std::setfill('0')
        << elapsed() << " sec" << std::endl;
    std::cout << oss.str();
}

double TicToc::elapsed() const {
    if (!is_working_)
        return dur_.count();
    else
        return (dur_ + std::chrono::duration_cast<duration_t>(hr_clock_t::now() - tp_)).count();
}
