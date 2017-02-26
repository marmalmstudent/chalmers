#include "Handin5.h"
#include <iostream>
#include <string.h>

int main(int argc, char *argv[], char *envp[])
{
    if (argc > 1)
    {
        if (strcmp(argv[1], "1") == 0)
        {
            std::cout << "Initializing task 1" << std::endl;
            Handin5 ha5;
            ha5.initTask1();
        }
        else if (strcmp(argv[1], "2") == 0)
        {
            std::cout << "Initializing task 2" << std::endl;
        }
        else if (strcmp(argv[1], "3") == 0)
        {
            std::cout << "Initializing task 3" << std::endl;
        }
        else if (strcmp(argv[1], "4") == 0)
        {
            std::cout << "Initializing task 4" << std::endl;
        }
        else if (strcmp(argv[1], "5") == 0)
        {
            std::cout << "Initializing task 5" << std::endl;
        }
        else if (strcmp(argv[1], "6") == 0)
        {
            std::cout << "Initializing task 6" << std::endl;
        }
        else
        {
            std::cerr << "Usage: /path/to/file <task_nbr>" << std::endl;
            return -1;
        }
    }
    else
    {
        std::cerr << "Usage: /path/to/file <task_nbr>" << std::endl;
        return -1;
    }
}
