/*
 * stringutil.h
 *
 *  Created on: 2022年05月23日
 *      Author: jwzhu
 */

#ifndef STRINGUTIL_H_
#define STRINGUTIL_H_
#include <algorithm> 
#include <cctype>
#include <locale>
namespace NSPutils {

static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}
 
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}
 
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}
 
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}
 
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

/**
 * Remove all leading and trailing spaces from the input. The result is a trimmed copy of the input.
 */  
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}


static float toValue_impl(const std::string &str, float) { 
    return stof(str); 
}

static double toValue_impl(const std::string &str, double) { 
    return stod(str); 
}

static long double toValue_impl(const std::string &str, long double) { 
    return stod(str); 
}

static int toValue_impl(const std::string &str, int) { 
    return stoi(str); 
}

static long toValue_impl(const std::string &str, long) { 
    return stol(str); 
}

static long long toValue_impl(const std::string &str, long long) { 
    return stoll(str); 
}

static unsigned long toValue_impl(const std::string &str, unsigned long) { 
    return stoul(str); 
}

static unsigned long long toValue_impl(const std::string &str, unsigned long long) { 
    return stoull(str); 
}

/**
 * Convert a string into numerical values.
 */ 
template<class T> 
T toValue(const std::string &str) { 
    return toValue_impl(str, T()); 
}


/**
 * Split a string by char delimiters.
 */ 
static void split(const std::string& s, std::vector<std::string>& tokens, const std::string& delimiters){

	std::string::size_type start = s.find_first_not_of(delimiters, 0);
	std::string::size_type pos = s.find_first_of(delimiters, start);
	while (pos != std::string::npos || start != std::string::npos){
		tokens.emplace_back(s.substr(start, pos - start));
		start = s.find_first_not_of(delimiters, pos);
		pos = s.find_first_of(delimiters, start);
	}
}

/**
 * Split a string by column widths.
 */
static void split(const std::string& s, std::vector<std::string>& tokens, const std::vector<int> & cols){
    int start = 0;
    for (auto &col : cols){
        tokens.emplace_back(NSPutils::trim_copy(s.substr(start,col)));
        start = start + col;
    }
}

}

#endif /* STRINGUTIL_H_ */
