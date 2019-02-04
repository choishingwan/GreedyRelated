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

    /*!
     * \brief Function to check if two double are equal from
     *        https://stackoverflow.com/a/4010279/1441789
     * \param a the first double
     * \param b the second double
     * \param error_factor level of error, should be of no concern to us at the
     *        moment
     * \return True if two double are equal
     */
    inline bool logically_equal(double a, double b, double error_factor = 1.0)
    {
        return ((a == b)
                || (std::abs(a - b) < std::abs(std::min(a, b))
                                          * std::numeric_limits<double>::epsilon()
                                          * error_factor));
    }
}
#endif /* misc_hpp */
