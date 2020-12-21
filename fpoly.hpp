#ifndef _FPOLY_HPP_
#define _FPOLY_HPP_

#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>

#include "globle.hpp"

#define M_PI 3.1415926535897932384

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

typedef std::complex<double> cp;

class fpolynome:public std::vector<double>
{
    private:
        //std::vector<double> v;
        
        vector<cp> fft(vector<cp> a,bool f);

    public:
        int n;
        fpolynome(int _n=cn): n(_n),vector<double>(2048) {}

        void ref(void);
        void inref(void);

        //double& operator[](int t) {return v[t];}

        friend std::ostream& operator<<(std::ostream& out,fpolynome a);
        
        friend fpolynome operator+(fpolynome a,fpolynome b);
        friend fpolynome operator-(fpolynome a,fpolynome b);
        friend fpolynome operator+(fpolynome a,double b);
        friend fpolynome operator-(fpolynome a,double b);
        friend fpolynome operator*(fpolynome a,double b);
        friend fpolynome operator*(double a,fpolynome b);        
        friend fpolynome operator*(fpolynome a,fpolynome b);   

        double ans(double x);
};

void fpolynome::ref(void)
{
    for(vector<double>::iterator i=begin()+n/2,j=begin();i<end();i++,j++)
    {
        *j-=*i;
        *i=0;
    }
}

void fpolynome::inref(void)
{
    ref();
}

std::vector<cp> fpolynome::fft(vector<cp> a,bool f)
{
    getrev(n);

    for(int i=0;i<n;i++)
    {
        if(i<nttrev[i])
        {
            std::swap(a[i],a[nttrev[i]]);
        }
    }
    //std::cout<<"Asdf";
    int ind=int(f)*2-1;
    for(int i=1;i<n;i*=2)
    {
        cp temp(cos(M_PI/i),sin(M_PI/i)*ind);
        for(int j=0;j<n;j+=i*2)
        {
            cp omega(1,0);
            for(int k=0;k<i;k++,omega*=temp)
            {
                // std::cout<<j+k<<" "<<j+k+i<<"$"<<std::endl;
                cp x=a[j+k],y=omega*a[j+k+i];

                a[j+k]=x+y;
                a[j+k+i]=x-y;
            }
        }
    }

    return a;
}

std::ostream& operator<<(std::ostream& out,fpolynome a)
{
    for(int i=0;i<a.n;i++)
    {
        out<<a[i]<<" ";
    }
    return out;
}
 
fpolynome operator+(fpolynome a,fpolynome b)
{
    fpolynome t(std::max(a.n,b.n));

    for(int i=0;i<a.n;i++) t[i]=a[i]+b[i];
    t.ref();

    return t;
}
 
fpolynome operator-(fpolynome a,fpolynome b)
{
    fpolynome t(std::max(a.n,b.n));

    for(int i=0;i<a.n;i++) t[i]=a[i]-b[i];
    t.ref();

    return t;
}

fpolynome operator+(fpolynome a,double b)
{
    fpolynome t=a;

    t[0]=a[0]+b;
    t.ref();

    return t;
}

fpolynome operator-(fpolynome a,double b)
{
    fpolynome t=a;

    t[0]=a[0]-b;
    t.ref();

    return t;
}
 
fpolynome operator*(fpolynome a,double b)
{
    fpolynome t(a.n);

    for(int i=0;i<a.n;i++) t[i]=a[i]*b;
    t.ref();

    return t;
}

fpolynome operator*(double a,fpolynome b)
{
    fpolynome t(b.n);

for(int i=0;i<b.n;i++) t[i]=a*b[i];

    t.ref();

    return t;
}
 
fpolynome operator*(fpolynome a,fpolynome b)
{
    fpolynome c;

    std::vector<cp> *ta=new std::vector<cp>,*tb=new std::vector<cp>,*tc=new std::vector<cp>;
    for(int i=0;i<cn;i++)
    {
        ta->push_back({a[i],0.0});
        tb->push_back({b[i],0.0});
    }
    *ta=a.fft(*ta,true);
    *tb=a.fft(*tb,true);

    for(int i=0;i<2048;i++)
    {
        tc->push_back((*ta)[i]*(*tb)[i]);
    }

    *tc=c.fft(*tc,false);
    for(int i=0;i<cn;i++)
    {
        c[i]=(*tc)[i].real();
    }

    return c;
}

double fpolynome::ans(double x)
{
    double tans=0;

    for(std::vector<double>::reverse_iterator i=rbegin();i<rend();i++)
    {
        tans*=x;
        tans+=*rbegin();
    }

    return tans;
}

fpolynome inv(fpolynome f,int n)
{
    int k=0,degf,degg=n;
    fpolynome b,c,g;
    b[0]=1.0;
    g[0]=1.0,g[n]=1.0;
    for(degf=cn-1;degf>=0;degf--)
    {
        if(f[degf]!=0.0) break;
    }

    while(true)
    {
        while(f[0]==0.0)
        {
            if(degf<=0) return b;
            for(int i=1;i<cn;i++)
            {
                f[i-1]=f[i];                
            }
            f[cn-1]=0.0;
            for(int i=cn-1;i>0;i--)
            {
                g[i]=g[i-1];
            }
            g[0]=0.0;
            k++;
            degf--;
            degg++;
        }
        if(degf==0)
        {
            b=b*(1.0/f[0]);
            for(int i=cn-1;i>=0;i--)
            {
                if(i>=k) b[i]=b[i-k];
                else b[i]=0.0;
            }
            b.n=n;
            b.ref();
            return b;
        }
        if(degf<degg)
        {
            swap(f,g);
            swap(b,c);
        }
        double u=f[0]/g[0];
        f=f-u*g;
        b=b-u*c;
    }
}

fpolynome trans(fpolynome a)
{
    for(int i=1;i<a.n;i++)
    {
        if(i<a.n-i) 
        {
            std::swap(a[i],a[a.n-i]);
            a[i]=-a[i];
            a[a.n-i]=-a[a.n-i];
        }
        else if(i==a.n-i) a[i]=-a[i];
    }

    return a;
}

#endif