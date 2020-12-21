#ifndef _RANDOM_HPP_
#define _RANDOM_HPP_

#include <random>
#include <ctime>
#include <cmath>

#include "globle.hpp"
#include "poly.hpp"

typedef long long ll;

double normrho[2048];

polynome rand_u(int n=cn,int c=0,int q=cq)
{
    std::random_device rd;
    std::default_random_engine e(rd());
    std::uniform_int_distribution<int> u(q/2-q,q/2);

    polynome ans(n,q);
    for(int i=0;i<n;i++)
    {
        ans[i]=rd()%q;
        if(ans[i]>q/2) ans[i]-=q;
    }

    return ans;
}
polynome rand_d(int n=cn,double c=0.0,double sigma=csigma)
{
    std::random_device rd;
    std::default_random_engine e(rd());
    std::normal_distribution<> nd(c,sigma);

    int maxa=int(round(nd(e)))*(int(rd())%2*2-1);

    std::default_random_engine e2(rd());
    std::uniform_int_distribution<int> u(-maxa,maxa);

    polynome ans(n);

    for(int i=0;i<n;i++)
    {
        ans[i]=u(e2);
    }
    ans[rd()%n]=maxa;

    return ans;
}


std::vector<double> rand_p(double sig,std::vector<double>& l,std::vector<double>& h)
{
    std::random_device rd;
    std::default_random_engine e(rd());
    std::normal_distribution<> nd(0.0,sig);

    std::vector<double> z,p;

    for(int i=0;i<l.size();i++)
    {
        z.push_back(nd(e));
    }
    for(int i=0;i<l.size()-1;i++)
    {
        p.push_back(l[i]*z[i]+h[i+1]*z[i+1]);
    }
    p.push_back(*h.rbegin()*(*z.rbegin()));

    return p;
}
polynome rand_d(double sig,std::vector<double> c,std::vector<double> d)
{
    std::random_device rd;
    std::default_random_engine e(rd());
    std::normal_distribution<> nd1(-(*c.rbegin())/(*d.rbegin()),sig/(*d.rbegin()));

    polynome z(c.size());

    *z.rbegin()=ll(round(nd1(e)));
    for(int i=0;i<c.size();i++)
    {
        c[i]-=*z.rbegin()*d[i];
    }
    for(int i=0;i<c.size()-1;i++)
    {
        std::normal_distribution<> nd(-c[i],sig);
        z[i]=ll(round(nd(e)));
    }

    return z;
}


#endif