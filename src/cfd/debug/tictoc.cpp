#include "tictoc.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace cfd::dbg;

namespace {

static constexpr size_t MIN_NAME_SIZE = 10;
static constexpr size_t REPORT_COLUMN_WIDTH = 15;

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
    virtual ~TicToc() = default;

    /// resets timer. Doesn't restarts it.
    void init();
    /// start timer
    void tic();
    /// pause timer
    void toc();
    /// print report to std::cout
    virtual void report() const;

    /// stop and report to std::cout
    void fintoc() {
        toc();
        report();
    }

    /// get elapsed time in seconds
    size_t elapsed() const;

    size_t calls_count() const {
        return calls_count_;
    }
    std::string name() const {
        return name_;
    }
    virtual size_t get_longest_name_size() const {
        return name_.size();
    }

private:
    using hr_clock_t = std::chrono::high_resolution_clock;
    using time_point_t = hr_clock_t::time_point;
    using duration_t = std::chrono::milliseconds;

    std::string name_;
    bool is_working_;
    time_point_t tp_;
    duration_t dur_;
    size_t calls_count_ = 0;
};

class IndexedTicToc : public TicToc {
public:
    /// create time with defined name.
    IndexedTicToc(std::string name) : TicToc(name, false){};

    void tic_indexed(size_t index) {
        tic();
        auto fnd = entries_.find(index);
        if (fnd == entries_.end()) {
            std::string new_name = "" + name() + "-" + std::to_string(index);
            entries_.emplace(index, TicToc(new_name, true));
        } else {
            fnd->second.tic();
        }
    }

    void toc_indexed(size_t index) {
        toc();
        auto fnd = entries_.find(index);
        if (fnd != entries_.end()) {
            fnd->second.toc();
        }
    }

    void report() const override {
        TicToc::report();
        for (auto& [_, v]: entries_) {
            v.report();
        }
    }

    size_t get_longest_name_size() const override {
        size_t ret = name().size();
        for (auto& [_, v]: entries_) {
            ret = std::max(ret, v.get_longest_name_size());
        }
        return ret;
    }

private:
    std::map<size_t, TicToc> entries_;
};

struct Timers {
    std::map<std::string, std::shared_ptr<TicToc>> data;
    size_t longest_name_size = MIN_NAME_SIZE;

    ~Timers() {
        FinReport();
    }

    bool has(std::string s) {
        return data.find(s) != data.end();
    }

    TicToc& get(std::string s) {
        auto fnd = data.find(s);
        if (fnd == data.end()) {
            auto eres = data.emplace(s, new TicToc(s, 0));
            return *eres.first->second;
        } else {
            return *fnd->second;
        }
    }

    IndexedTicToc& get_indexed(std::string s) {
        auto fnd = data.find(s);
        if (fnd == data.end()) {
            auto eres = data.emplace(s, new IndexedTicToc(s));
            return *dynamic_cast<IndexedTicToc*>(eres.first->second.get());
        } else {
            auto ptr = dynamic_cast<IndexedTicToc*>(fnd->second.get());
            if (ptr != nullptr) {
                return *ptr;
            } else {
                throw std::runtime_error("Timer names " + s + " is not an indexed timer");
            }
        }
    }

    std::vector<std::string> keys() {
        std::vector<std::string> ret;
        for (auto& v: data) {
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

    void recalculate_longest_name() {
        longest_name_size = MIN_NAME_SIZE;
        for (const auto& [_, v]: data) {
            longest_name_size = std::max(longest_name_size, v->get_longest_name_size());
        }
    }

    size_t size() const {
        return data.size();
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

void cfd::dbg::Tic(std::string s, size_t index) {
    auto& tm = alltimers_.get_indexed(s);
    tm.tic_indexed(index);
}

void cfd::dbg::Tic1(std::string s) {
    cfd::dbg::Toc();
    cfd::dbg::Tic(s);
}

void cfd::dbg::Toc(std::string s) {
    if (s.size() == 0) {
        for (auto k: alltimers_.keys())
            Toc(k);
    } else {
        if (alltimers_.has(s)) {
            alltimers_.get(s).toc();
        };
    }
}

void cfd::dbg::Toc(std::string s, size_t index) {
    auto& tm = alltimers_.get_indexed(s);
    tm.toc_indexed(index);
}

void cfd::dbg::Report(std::string s) {
    if (s.size() == 0) {
        for (auto k: alltimers_.keys())
            Report(k);
    } else {
        if (alltimers_.has(s)) {
            alltimers_.get(s).report();
        }
    }
}

void cfd::dbg::FinReport(std::string s) {
    if (s.size() == 0) {
        if (alltimers_.size() == 0) {
            return;
        }
        Toc();
        alltimers_.recalculate_longest_name();
        std::ostringstream oss;
        oss << "------ Tictoc report (ms)" << std::endl;
        oss << std::left << std::setw(static_cast<int>(alltimers_.longest_name_size)) << "Process";
        oss << std::right << std::setw(REPORT_COLUMN_WIDTH) << "Total time";
        oss << std::right << std::setw(REPORT_COLUMN_WIDTH) << "Calls count";
        oss << std::right << std::setw(REPORT_COLUMN_WIDTH) << "Time per call";
        oss << std::endl;
        for (size_t _ = 0; _ < alltimers_.longest_name_size + 3 * REPORT_COLUMN_WIDTH; ++_) {
            oss << '-';
        }
        oss << std::endl;
        std::cout << oss.str();
        for (auto k: alltimers_.keys_sorted_by_elapsed()) {
            FinReport(k);
        }
    } else {
        if (alltimers_.has(s)) {
            alltimers_.get(s).fintoc();
        }
        alltimers_.erase(s);
    }
}

TicToc::TicToc(std::string name, bool start) : name_(name), is_working_(false), dur_(duration_t::zero()) {
    if (start)
        tic();
}

void TicToc::tic() {
    if (!is_working_) {
        is_working_ = true;
        calls_count_ += 1;
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
    if (calls_count_ == 0) {
        return;
    }
    std::ostringstream oss;
    oss << std::left << std::setw(static_cast<uint32_t>(alltimers_.longest_name_size)) << name_ << std::right
        << std::setw(REPORT_COLUMN_WIDTH) << elapsed() << std::right << std::setw(REPORT_COLUMN_WIDTH) << calls_count()
        << std::right << std::setw(REPORT_COLUMN_WIDTH) << elapsed() / calls_count();
    oss << std::endl;
    std::cout << oss.str();
}

size_t TicToc::elapsed() const {
    if (!is_working_)
        return dur_.count();
    else
        return (dur_ + std::chrono::duration_cast<duration_t>(hr_clock_t::now() - tp_)).count();
}

ScopedTic::ScopedTic(std::string s) : s_(s), is_indexed_(false) {
    Tic(s);
}
ScopedTic::ScopedTic(std::string s, size_t index) : s_(s), is_indexed_(true), index_(index) {
    Tic(s, index);
}
ScopedTic::~ScopedTic() {
    if (!is_indexed_) {
        Toc(s_);
    } else {
        Toc(s_, index_);
    }
}
