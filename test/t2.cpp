#include <iostream>
#include "Aggrate/Memory.h"

int main()
{
    MemoryArena area;
    int *p = area.Alloc<int>(100);
        

    return 0;
}