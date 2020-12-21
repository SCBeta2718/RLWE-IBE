#ifndef _CLIENT_HPP_
#define _CLIENT_HPP_

#include <string>

#include "trapdoor.hpp"
#include "poly.hpp"

class client
{
    private:
        trapdoor sv;
        polynome beta;
        ll q;
    
    public:
        polynome s;
        matrix e0;
        matrix c0;
        polynome e1;
        polynome c1;
        matrix omega;

        client(std::string id,double sigma): sv(),q(cq),beta(id_hash(id,cq)),s(rand_u(cn/2,0,q)),e1(rand_d(cn/2,0.0,sigma))
        {
            for(int i=0;i<sv.A.size();i++)
            {
                e0.push_back(rand_d(cn/2,0.0,sigma));
            }
            for(int i=0;i<sv.A.size();i++)
            {
                c0.push_back(e0[i]+sv.A[i]*s);
            }
            c1=beta*s+e1;

            omega=gausssamp(sv.A,sv.rho,sv.v,beta,csigma);
        }
        polynome enc(std::string s)
        {
            polynome miu(cn,cq);
            for(int i=0;i<s.length();i++)
            {
                miu[i]=(s[i]-'0')*((q+1)/2);
            }

            return c1+miu;
        }
};

std::string dec(matrix c0,polynome c1,matrix omega)
{
    std::string ans;
    for(int i=0;i<omega.size();i++)
    {
        c1=c1-omega[i]*c0[i];
    }
    for(int i=0;i<cn;i++)
    {
        if(abs(c1[i])<cq/4) ans+='0';
        else ans+='1';
    }

    return ans;
}

#endif
