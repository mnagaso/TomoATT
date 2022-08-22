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

    void start_timer(){
        start = std::clock();
    };
    void stop_timer(){
        end = std::clock();

        // print time
        if (if_verbose)
            std::cout << "Time for " << name << ": " << (double)(end - start) / CLOCKS_PER_SEC << " sec" << std::endl;
    };
};


#endif // TIMER_H