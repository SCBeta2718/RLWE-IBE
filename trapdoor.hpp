#ifndef _SERVER_HPP_
#define _SERVER_HPP_

#include <vector>
#include <algorithm>

#include "globle.hpp"
#include "poly.hpp"
#include "fpoly.hpp"
#include "random.hpp"
#include "sampling.hpp"

class trapdoor
{
    private:
        int sigma,q,n,_m;
 
    public:
        matrix A;
        matrix rho,v;

        trapdoor(int _n=cn/2,int _q=cq,double _sigma=csigma): n(_n),q(_q),sigma(_sigma),_m(ilog2(_q)+1)
        {
            polynome a(rand_u(n,q));
            
            for(int i=0;i<_m;i++)
            {
                rho.push_back(rand_d(n,0.0,sigma));
            }
            for(int i=0;i<_m;i++)
            {
                v.push_back(rand_d(n,0.0,sigma));
            }

            A.push_back(a);
            A.push_back(polynome(n,q));
            A[1][0]=1;

            for(int i=0;i<_m;i++)
            {
                A.push_back(polynome(n,q)-(a*rho[i]+v[i])+qpow(2,i,q));
            }
        }
};


#endif

