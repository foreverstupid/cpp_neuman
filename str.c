#include "str.h"

/* convert string to integer */
int str2int(const char *str)
{
    int i;
    int res = 0;

    if(!str){
        return 0;
    }

    for(i = 0; str[i]; i++){
        res *= 10;
        res += str[i] - '0';
    }

    return res;
}



/* convert string to double */
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
