#include "string_operations.hpp"

bool isNumber(const char *str)
{
    int i = 0;
    bool wasPoint = false;

    if(str[0] == '-'){
        i = 1;
    }

    for(; str[i]; i++){
        if(str[i] == '.'){
            if(!wasPoint){
                wasPoint = true;
                continue;
            }else{
                return false;
            }
        }

        if(str[i] < '0' || str[i] > '9'){
            return false;
        }
    }

    return true;
}



int str2int(const char *str)
{
    int res = 0;
    int i = 0;
    int sign = 1;

    if(!str){
        return 0;
    }

    if(str[0] == '-'){
        sign = -1;
        i = 1;
    }

    for(; str[i]; i++){
        res *= 10;
        res += str[i] - '0';
    }

    return res * sign;
}



double str2double(const char *str)
{
    int i = 0;
    int sign = 1;
    double res = 0.0;
    double hlp = 0.1;

    if(!str){
        return 0.0;
    }

    if(str[0] == '-'){
        sign = -1;
        i = 1;
    }

    for(; str[i] && str[i] != '.'; i++){
        res *= 10;
        res += str[i] - '0';
    }

    if(str[i] == '.'){
        for(++i; str[i]; i++){
            res += hlp * (str[i] - '0');
            hlp /= 10;
        }
    }

    return res * sign;
}
