//
//  misc.hpp
//  plink
//
//  Created by Shing Wan Choi on 18/08/2016.
//  Copyright © 2016 Shing Wan Choi. All rights reserved.
//

#ifndef misc_hpp
#define misc_hpp

#include <stdio.h>
#include <stdexcept>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

namespace misc{
    
    //Functions from R
    double dnorm(double x, double mu=0.0, double sigma=1.0, bool log=false);
    double qnorm(double p, double mu=0.0, double sigma=1.0, bool lower_tail=true, bool log_p=false);
    
    // codes from stackoverflow
    std::vector<std::string> split(const std::string seq, const std::string separators="\t ");
    template <typename T> inline
    T convert(const std::string& str){
        std::istringstream iss(str);
        T obj;
        
        iss >> std::ws >> obj >> std::ws;
        
        if(!iss.eof())
            throw std::runtime_error("Unable to convert the input");
        
        return obj;
    }
    // trim from start (in place)
    void ltrim(std::string &s);
    // trim from end (in place)
    void rtrim(std::string &s);
    // trim from both ends (in place)
    void trim(std::string &s);
    // trim from start (copying)
    std::string ltrimmed(std::string s);
    // trim from end (copying)
    std::string rtrimmed(std::string s);
    // trim from both ends (copying)
    std::string trimmed(std::string s);

}
#endif /* misc_hpp */
