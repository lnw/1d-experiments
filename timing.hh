#ifndef TIMING_HH
#define TIMING_HH

#include <vector>
#include <ostream>
#include <fstream>
#include <chrono>

using namespace std;
typedef std::chrono::duration<long int, ratio<1, 1000000000>> nano_s;
typedef std::chrono::high_resolution_clock Clock;

class time_item{
  string label;
  unsigned indent;
  nano_s duration; // in ms

public:
  time_item(const string& l, unsigned i, nano_s d): label(l), indent(i), duration(d) {};

  friend ostream& operator<<(ostream& S, const time_item& ti) {
    for(int i=0; i<ti.indent; i++) S << "  ";
    S << ti.label << ": ";
    S << std::chrono::duration_cast<std::chrono::milliseconds>(ti.duration).count() << " ms";
    return S;
  }
};


class timing{
  vector<time_item> times;
  unsigned indentation_state;

public:

  timing(): indentation_state(0) {};

  ~timing(){
    ofstream log;
    log.open ("timings.log", ios::ate);
    for(int i=0; i<times.size(); i++)
      log << times[i] << endl;
    log.close();
  }

   void ind_up(){indentation_state++;}
   void ind_down(){indentation_state--;}

   void append_item(const char* l, nano_s d){
     times.push_back(time_item(string(l), indentation_state, d));
   }
};

#endif

