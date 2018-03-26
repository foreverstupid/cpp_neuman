#ifndef STRING_OPERATIONS_MODULE_H
#define STRING_OPERATIONS_MODULE_H

#ifdef __cplusplus
extern "C"
{
#endif

int str2int(const char *str);

double str2double(const char *str);

static inline int is_num_arg(const char *str)
{
    return str && ((str[0] >= '0' && str[0] <= '9') || str[0] == '.');
}

#ifdef __cplusplus
}
#endif

#endif
