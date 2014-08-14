
/*test.cc*/
#include <iostream>
#include "plotGrid.h"

int main() 
{

    int nx = 200;
    int* states = new int[nx*nx];

    for (int i=0;i<nx;i++)
        for (int j=0;j<nx;j++)
        {
            states[j*nx+i] = (i%2)*(j%2);
        }

    plotGrid* pg = new plotGrid;

    for (int t=0;t<100;t++)
        pg->draw(nx, states);
    return 0;
}
