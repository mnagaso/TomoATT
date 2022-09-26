#ifndef TIMER_H
#define TIMER_H

#include<iostream>
#include<ctime>
#include<chrono>
#include<string>

class Timer {
public:
    Timer(std::string& name){
        this->name = name;
        this->start_timer();
    };
    ~Timer(){};

    std::string name;
    std::clock_t start;
    std::clock_t end;
    std::clock_t tmp_delta;
    std::clock_t time_back;

    void start_timer(){
        start = std::clock();
        time_back = start;
    };
    void stop_timer(){
        end = std::clock();

        // print time
        if (if_verbose)
            std::cout << "Time for " << name << ": " << (double)(end - start) / CLOCKS_PER_SEC << " sec" << std::endl;
    };

    double get_t_delta(){
        tmp_delta = std::clock() - time_back;
        time_back = std::clock();
        return (double) tmp_delta / CLOCKS_PER_SEC;
    };
};


#endif // TIMER_H