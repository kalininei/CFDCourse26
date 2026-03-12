#include "cfd/debug/printer.hpp"
#include "cfd/debug/tictoc.hpp"
#include "cfd/mat/lodmat.hpp"
#include "cfd26_test.hpp"
#include <chrono>
#include <thread>

using namespace cfd;

TEST_CASE("Matrix printer test", "[matrix-printer]") {
    LodMatrix M(3);
    M.add_value(0, 0, 1.0);
    M.add_value(1, 0, 2.0);
    M.add_value(1, 1, 4.0);
    M.add_value(2, 2, 3.0);
    dbg::print(M);
    dbg::print(1, M);
}

TEST_CASE("Tictoc test", "[tictoc]") {
    auto tglob = dbg::ScopedTic("ALL");
    {
        auto t = dbg::ScopedTic("MYPROC12", 121);
        dbg::Tic("process1");
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        dbg::Toc("process1");

        dbg::Tic("process2");
        std::this_thread::sleep_for(std::chrono::milliseconds(6));
        dbg::Toc("process2");
    }

    dbg::Tic("process1");
    std::this_thread::sleep_for(std::chrono::milliseconds(4));
    dbg::Toc("process1");

    dbg::Tic("process1");
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    dbg::Toc("process1");

    for (size_t j = 1; j < 2; ++j) {
        for (size_t i = 0; i < 3; ++i) {
            dbg::Tic("loop", i);
            std::this_thread::sleep_for(std::chrono::milliseconds(5 + i));
            dbg::Toc("loop", i);
        }
    }

    {
        auto tic = dbg::ScopedTic("process", 1);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    {
        auto tic = dbg::ScopedTic("process", 2);
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    {
        auto tic = dbg::ScopedTic("process", 3);
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    dbg::FinReport();
}
