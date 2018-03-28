#ifndef STRING_OPERATIONS_MODULE_HPP
#define STRING_OPERATIONS_MODULE_HPP

inline bool isDigit(char ch)
{
    return ch >= '0' && ch <= '9';
}

bool isNumber(const char *str);

int str2int(const char *str);

double str2double(const char *str);

#endif
