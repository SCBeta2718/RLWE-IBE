#include <iostream>
#include <cstdlib>

#include "client.hpp"

using namespace std;


int main(void)
{
    system("ulimit -s 524288");

    polynome a=rand_d(cn/2,0.0,10*csigma);

    string id,miu0,miu1;
    cout<<"input id:\n";
    cin>>id;
    client alice(id,csigma);
    cout<<"input 01 string:\n";
    cin>>miu0;
    polynome c1=alice.enc(miu0);
    cout<<"============\n";
    miu1=dec(alice.c0,c1,alice.omega);
    cout<<miu1;
    
    return 0;
}