#include <chrono>
#include <iostream>
#include <cmath>
#include <iomanip>

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
    
    void update_progress(const float progress, const string extra_string = ""){
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

            cout << "\r\033[F" << "ETA: " << ETA_h << "h  " << ETA_m << "m  " << ETA_s << "s  \t\t" << extra_string << "\n" << "[" << string(pos,'=') << string(bar_width-pos,' ') << "]" << setprecision(3)<< (progress * 100.0) << "\t%" << flush;
        }else
        {
            cout << "Invalid Value provided progress cannot exceed 100%\r" << flush;

        }
    }
    
    
    ~ProgressBar(){};
};
