#include <chrono>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include <stdio.h>

class ProgressBar
{
private:
    typedef std::chrono::system_clock::time_point time_t;
    std::chrono::system_clock clock = std::chrono::system_clock();
    char progress_symbol = '#';
    int bar_length = 50;
    time_t t_start = clock.now();
    time_t t1   = clock.now();
    time_t t2   = clock.now();
    float p1 = 0.f;
    float p2 = 0.f;
    
public:
    ProgressBar(){};
    ProgressBar(char progress_symbol, int bar_length): progress_symbol(progress_symbol), bar_length(bar_length){}
    
    void update_progress(const float progress, const std::string extra_string = "", std::vector<std::pair<std::string, std::chrono::nanoseconds>> extra_timers = {}){
        using namespace std::chrono;

        if(progress <= 1){
            p1 = p2;
            p2 = progress;

            t1 = t2;
            t2 = clock.now();

            int ETA = floor((1.f-progress) / ( (progress) / ((t2-t_start).count()/1e9)) );
            int ETA_h = ETA / 3600;
            int ETA_m = (ETA - (ETA_h * 3600))/60;
            int ETA_s = ETA - (ETA_h * 3600) - (ETA_m*60);

            int bar_width = 70;
            int pos = bar_width * progress;
            //q: what does \r\033[F do?
            //a: \r moves the cursor to the beginning of the line, \033[F moves the cursor up one line
            
            
            std::cout << "\r\033[F" << "ETA: " << std::setw(3) << ETA_h << "h  "<< std::setw(3) << ETA_m << "m  " << std::setw(3) << ETA_s << "s  " << std::setw(30) << extra_string << "\n" << "[" << std::string(pos,'=') << std::string(bar_width-pos,' ') << "]" << std::setprecision(3)<< std::setw(7) << (progress * 100.0) 
            << "%    \n";
            for(auto timer : extra_timers){
                std::cout << std::setw(10) << timer.first << ": " << std::setprecision(2) <<  std::setw(5) << (100.f*float(timer.second.count()) / float((t2-t_start).count())) << "%   \n";
            }
            for (size_t i = 0; i < extra_timers.size() + 1; i++)
            {
                std::cout<< "\033[F";
            }
            std::cout << std::flush;
        }else
        {
            std::cout << "Invalid Value provided progress cannot exceed 100%\r" << std::flush;

        }
    }
    
    
        ~ProgressBar() {
            std::cout << "\n" << "\r\033[K" << std::flush; //clear line
        }
    };
