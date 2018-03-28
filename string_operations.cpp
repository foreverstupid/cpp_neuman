#include "string_operations.hpp"

bool isNumber(const char *str)
{
    bool wasPoint = false;

    for(int i = 0; str[i]; i++){
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

    if(!str){
        return 0;
    }

    for(int i = 0; str[i]; i++){
        res *= 10;
        res += str[i] - '0';
    }

    return res;
}



double str2double(const char *str)
{
    int i;
    double res = 0.0;
    double hlp = 0.1;

    if(!str){
        return 0.0;
    }

    for(i = 0; str[i] && str[i] != '.'; i++){
        res *= 10;
        res += str[i] - '0';
    }

    if(str[i] == '.'){
        for(++i; str[i]; i++){
            res += hlp * (str[i] - '0');
            hlp /= 10;
        }
    }

    return res;

}
