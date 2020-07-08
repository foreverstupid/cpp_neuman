#ifndef STRING_OPERATIONS_MODULE_HPP
#define STRING_OPERATIONS_MODULE_HPP

#include <stdio.h>

inline bool isDigit(char ch)
{
    return ch >= '0' && ch <= '9';
}

bool isNumber(const char *str);

int str2int(const char *str);

double str2double(const char *str);

bool equals(const char *str1, const char *str2);

#endif
