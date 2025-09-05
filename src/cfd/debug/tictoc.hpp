/**
 * @file
 *
 * @brief timer for quick profiling
 */
#ifndef CFD_TICTOC_HPP
#define CFD_TICTOC_HPP

#include <chrono>
#include <string>

namespace cfd{
namespace dbg{

/**
 * @brief starts a global timer with string id.
 * @param s string identifier of this timer. Default is "timer1,2,3, etc"
 */
void Tic(std::string s = "");

/**
* @brief starts a global timer and stops all others
*/
void Tic1(std::string s = "");

/**
 * @brief stops global timers
 * @param s timer idetifier. By default stops all timers
 */
void Toc(std::string s = "");

/**
 * @brief prints global timer report to std::cout
 * @param s timer identifier. By default reports all defined timers
 */
void Report(std::string s = "");

/**
 * @brief reports and deletes global timer
 * @param s timer identifer. By default reports and removes all timers
 *
 * This procedure is always called by program termination.
 * So there is no need to call it explicitly.
 */
void FinReport(std::string s = "");

} // namespace dbg
} // namespace cfd

#endif
