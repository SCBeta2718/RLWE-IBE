#ifndef _POLY_HPP_
#define _POLY_HPP_

#include <vector>
#include <iostream>
#include <algorithm>

#include "globle.hpp"


#ifndef _REV_
#define _REV_

int nttrev[cn];

void getrev(int n=cn)
{
    int l=ilog2(n);
    for(int i=0;i<n;i++)
    {
        nttrev[i]=(nttrev[i/2]/2)|((i%2)<<(l-1));
    }
}

#endif

class polynome:public std::vector<ll>
{
    private:
        //std::vector<ll> v;
        
        polynome ntt(bool f);

    public:
        int n;
        ll q;
        polynome(int _n=cn,ll _q=-1): n(_n),q(_q),vector<ll>(2048) {}

        void ref(void);
        void inref(void);

        //ll& operator[](int t) {return v[t];}

        polynome ntt(void);
        polynome intt(void);

        friend std::ostream& operator<<(std::ostream& out,polynome a);
        
        friend polynome operator+(const polynome& a,const polynome& b);
        friend polynome operator-(const polynome& a,const polynome& b);
        friend polynome operator+(const polynome& a,ll b);
        friend polynome operator-(const polynome& a,ll b);
        friend polynome operator*(const polynome& a,ll b);
        friend polynome operator*(ll a,const polynome& b);        
        friend polynome operator*(const polynome& a,const polynome& b);
        friend polynome operator%(const polynome& a,ll b);       

        ll ans(ll x);
};

void polynome::ref(void)
{
    if(q==-1) return ;

    for(vector<ll>::iterator i=begin()+n/2,j=begin();i<end();i++,j++)
    {
        *j-=*i;
        *i=0;
    }
    for(vector<ll>::iterator i=begin();i<end();i++)
    {
        while(*i>q/2)
        {
            *i-=q;
        }
        while(*i<q/2-q)
        {
            *i+=q;
        }
    }
}

void polynome::inref(void)
{
    if(q==-1) return ;
    ref();
    for(std::vector<ll>::iterator i=begin();i<end();i++) *i=((*i)+q)%q;
}

polynome polynome::ntt(bool f)
{
    polynome a(*this);

    getrev(n);

    for(int i=0;i<n;i++)
    {
        if(i<nttrev[i])
        {
            std::swap(a[i],a[nttrev[i]]);
        }
    }
    for(int i=1;i<n;i*=2)
    {
        ll gn=qpow(cg,(q-1)/(i*2),q);
        for(int j=0;j<n;j+=i*2)
        {
            ll tg=1;
            for(int k=0;k<i;k++,tg=tg*gn%q)
            {
                // std::cout<<j+k<<" "<<j+k+i<<"$"<<std::endl;
                ll x=a[j+k],y=tg*a[j+k+i]%q;

                a[j+k]=(x+y)%q;
                a[j+k+i]=(x-y+q)%q;
            }
        }
    }

    if(f) return a;

    ll nv=qpow(n,q-2,q); 
    reverse(a.begin()+1,a.begin()+n);
    for(int i=0;i<n;i++)
    {
        a[i]=a[i]*nv%q;
    }

    return a;
}

polynome polynome::ntt(void) {return ntt(true);}

polynome polynome::intt(void) {return ntt(false);}

std::ostream& operator<<(std::ostream& out,polynome a)
{
    for(int i=0;i<a.n;i++)
    {
        out<<a[i]<<" ";
    }
    return out;
}
 
polynome operator+(const polynome& a,const polynome& b)
{
    polynome t(std::max(a.n,b.n),a.q);

    for(int i=0;i<a.n;i++) t[i]=a[i]+b[i];
    t.ref();

    return t;
}
 
polynome operator-(const polynome& a,const polynome& b)
{
    polynome t(std::max(a.n,b.n),a.q);

    for(int i=0;i<a.n;i++) t[i]=a[i]-b[i];
    t.ref();

    return t;
}

polynome operator+(const polynome& a,ll b)
{
    polynome t=a;

    t[0]=a[0]+b;
    t.ref();

    return t;
}

polynome operator-(const polynome& a,ll b)
{
    polynome t=a;

    t[0]=a[0]-b;
    t.ref();

    return t;
}
 
polynome operator*(const polynome& a,ll b)
{
    polynome t(a.n,a.q);

    if(a.q==-1) for(int i=0;i<a.n;i++) t[i]=a[i]*b;
    else for(int i=0;i<a.n;i++) t[i]=((a[i]+a.q)%a.q*b)%a.q;
    t.ref();

    return t;
}

polynome operator*(ll a,const polynome& b)
{
    polynome t(b.n,b.q);

    if(b.q==-1) for(int i=0;i<b.n;i++) t[i]=a*b[i];
    else for(int i=0;i<b.n;i++) t[i]=((a+b.q)%b.q*b[i])%b.q;

    t.ref();

    return t;
}
 
polynome operator*(const polynome& a,const polynome& b)
{
    polynome ta=a,tb=b,tc;

    ta.inref();
    tb.inref();

    ta=ta.ntt();
    tb=tb.ntt();

    for(int i=0;i<2048;i++)
    {
        tc[i]=ta[i]*tb[i]%cq;
    }

    tc=tc.intt();
    tc.ref();

    return tc;
}

polynome operator%(const polynome& a,ll b)
{
    polynome ans(a.n,b);

    for(int i=0;i<a.n;i++)
    {
        ans[i]=a[i]%b;
    }
    ans.ref();

    return ans;
}

ll polynome::ans(ll x)
{
    ll tans=0LL;

    for(std::vector<ll>::reverse_iterator i=rbegin();i<rend();i++)
    {
        tans*=x;
        tans+=*rbegin();

        if(q>0) tans%=q;
    }

    return tans;
}

polynome id_hash(std::string str,ll q)
{
    polynome ans(cn/2,q);

    for(int i=0;i<str.size();i++)    
    {
        ans[i]=ll(str[i]);
    }

    ans.ref();

    return ans;
}



#endif