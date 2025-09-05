#include "tictoc.hpp"
#include <iostream>
#include <map>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <sstream>

using namespace cfd::dbg;

namespace{

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
class TicToc{
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
	void fintoc(){
		toc();
		report();
	}

	/// get elapsed time in seconds
	double elapsed() const;

private:
	using hr_clock_t = std::chrono::high_resolution_clock;
	using time_point_t = hr_clock_t::time_point;
	using duration_t = std::chrono::duration<double>;

	std::string _name;
	bool _is_working;
	time_point_t _tp;
	duration_t _dur;
};


struct Timers{
	std::map<std::string, TicToc> data;

	~Timers(){
		FinReport();
	}

	bool has(std::string s){
		return data.find(s) != data.end();
	}

	TicToc& get(std::string s){
		auto fnd = data.find(s);
		if (fnd == data.end()){
			data[s] = TicToc(s, 0);
		}
		return data[s];
	}

	std::vector<std::string> keys(){
		std::vector<std::string> ret;
		for (auto& v: data){
			ret.push_back(v.first);
		}
		return ret;
	}

	std::vector<std::string> keys_sorted_by_elapsed(){
		std::vector<std::string> ret = keys();
		std::sort(ret.begin(), ret.end(),
			[this](const std::string& key1, const std::string& key2)->bool{
				return get(key1).elapsed() > get(key2).elapsed();
			});
		return ret;
	}

	void erase(std::string s){
		auto fnd = data.find(s);
		if (fnd != data.end()) data.erase(fnd);
	}
};

Timers _alltimers;

} // namespace

void cfd::dbg::Tic(std::string s){
	if (s.size() == 0){
		for (int i = 0; i < 99999; ++i){
			std::string nm = "Timer" + std::to_string(i);
			if (!_alltimers.has(nm)) return Tic(nm);
		}
	} else {
		auto& tm = _alltimers.get(s);
		tm.tic();
	}
}

void cfd::dbg::Tic1(std::string s){
	cfd::dbg::Toc();
	cfd::dbg::Tic(s);
}

void cfd::dbg::Toc(std::string s){
	if (s.size() == 0){
		for (auto k: _alltimers.keys())
			Toc(k);
	} else {
		if (_alltimers.has(s)){
			_alltimers.get(s).toc();
		}
	}
}

void cfd::dbg::Report(std::string s){
	if (s.size() == 0){
		for (auto k: _alltimers.keys())
			Report(k);
	} else {
		if (_alltimers.has(s)){
			_alltimers.get(s).report();
		}
	}
}

void cfd::dbg::FinReport(std::string s){
	if (s.size() == 0){
		Toc();
		for (auto k: _alltimers.keys_sorted_by_elapsed()){
			FinReport(k);
		}
	} else {
		if (_alltimers.has(s)){
			_alltimers.get(s).fintoc();
		}
		_alltimers.erase(s);
	}
}

TicToc::TicToc(std::string name, bool start)
: _name(name),
  _is_working(start),
  _dur(duration_t::zero()){
	if (start) tic();
}

void TicToc::init(){
	_dur = duration_t::zero();
	_tp = hr_clock_t::now();
}

void TicToc::tic(){
	if (!_is_working){
		_is_working = true;
		_tp = hr_clock_t::now();
	}
}

void TicToc::toc(){
	if (_is_working){
		_is_working = false;
		_dur += std::chrono::duration_cast<duration_t>(hr_clock_t::now() - _tp);
	}
}

void TicToc::report() const{
	std::ostringstream oss;
	oss << std::setw(10) << _name << ":  "
	    << std::setprecision(3) << std::setw(5) << std::left << std::setfill('0')
	    << elapsed() << " sec" << std::endl;
	std::cout << oss.str();
}

double TicToc::elapsed() const{
	if (!_is_working) return _dur.count();
	else return (_dur + std::chrono::duration_cast<duration_t>(hr_clock_t::now() - _tp)).count();
}
