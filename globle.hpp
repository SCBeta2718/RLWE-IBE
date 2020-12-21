#ifndef _GLOBLE_HPP
#define _GLOBLE_HPP

typedef long long ll;

const int cn=2048;
const int cl=11;
const ll cq=1004535809LL;
const double cmiu=149.3;
const double cdelta=1.003941;
const double csigma=4.578;
const ll cg=3LL;

ll qpow(ll a,ll n,ll p)
{
    ll ans=1LL;

    for(;n;n/=2,a=a*a%p)
    {
        if(n%2)
        {
            ans=ans*a%p;
        }
    }

    return ans;
}
ll ilog2(ll n)
{
    int i=0;

    while(n>1)
    {
        n/=2;
        i++;
    }

    return i;
}


#endif