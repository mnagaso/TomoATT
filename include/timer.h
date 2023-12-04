#ifndef TIMER_H
#define TIMER_H


#include <iostream>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <string>

class Timer {
public:
    Timer(std::string name, bool show_start_time = false) : name(name), show_start_time(show_start_time) {
        start_timer();
    }

    ~Timer() {}

    std::string name;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_time;
    bool show_start_time = false;

    void start_timer() {
        start = std::chrono::high_resolution_clock::now();
        time_back = start;
        if (show_start_time) {
            std::cout << this->name << " started at " << get_current_utc_time() << std::endl;
        }
    }

    void stop_timer() {
        end = std::chrono::high_resolution_clock::now();
        elapsed_time = end - start;
        print_elapsed_time();
    }

    double get_start() {
        return std::chrono::duration<double>(start.time_since_epoch()).count();
    }

    double get_t_delta() {
        auto current_time = std::chrono::high_resolution_clock::now();
        auto delta = current_time - time_back;
        time_back = current_time;
        return std::chrono::duration<double>(delta).count();
    }

    double get_t() {
        return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
    }

    std::string get_start_t() {
        return get_current_utc_time(start);
    }

    std::string get_end_t() {
        return get_current_utc_time(end);
    }

    std::string get_elapsed_t() {
        return std::to_string(elapsed_time.count()) + " sec";
    }

    std::string get_utc_from_time_t(std::time_t time_t) {
        std::tm utc_tm = *std::gmtime(&time_t);

        std::ostringstream oss;
        oss << std::put_time(&utc_tm, "%F %T") << " UTC";
        return oss.str();
    }

private:
    std::chrono::high_resolution_clock::time_point time_back;

    void print_elapsed_time() {
        if (if_verbose)
            std::cout << "Time for " << name << ": " << elapsed_time.count() << " sec" << std::endl;
    }

    std::string get_current_utc_time(std::chrono::high_resolution_clock::time_point time = std::chrono::high_resolution_clock::now()) {
        std::time_t time_t = std::chrono::system_clock::to_time_t(time);
        std::tm utc_tm = *std::gmtime(&time_t);

        std::ostringstream oss;
        oss << std::put_time(&utc_tm, "%F %T") << " UTC";
        return oss.str();
    }

};



#endif // TIMER_H