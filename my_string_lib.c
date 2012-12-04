

#include <stdlib.h>
#include "my_string_lib.h"




int str_starts_with(char * str, char first) {
    return *str == first;
}

double str_to_double(char * str) {
    return atof(str);
}

int str_to_int(char * str) {
    return atoi(str);
}

int count_chars(char * str, char c) {
    int count = 0;
    while (*str) { 
        if (*str == c) {
            count++;
        }
        str++;
    }
    return count;
}




