#include "string_operations.hpp"

bool isNumber(const char *str)
{
    double tmp;
    return sscanf(str, "%lf", &tmp) == 1;
}



int str2int(const char *str)
{
    int res;
    sscanf(str, "%d", &res);
    return res;
}



double str2double(const char *str)
{
    double res;
    sscanf(str, "%lf", &res);
    return res;
}



bool equals(const char *str1, const char *str2)
{
    int i;
    for(i = 0; str1[i]; i++){
        if(str1[i] != str2[i]){
            return false;
        }
    }

    return str2[i] == 0;
}
