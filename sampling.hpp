#ifndef _SAMPLING_HPP_
#define _SAMPLING_HPP_

#include <vector>

#include "poly.hpp"
#include "fpoly.hpp"
#include "random.hpp"

typedef std::vector<polynome> matrix;

std::pair<polynome,polynome> sample2z(fpolynome,fpolynome,fpolynome,fpolynome,fpolynome);
polynome samplefz(fpolynome,fpolynome);

std::pair<polynome,polynome> sample2z(fpolynome a,fpolynome b,fpolynome d,fpolynome c0,fpolynome c1)
{
    polynome q1=samplefz(d,c1);
    fpolynome _q1(q1.n);

    for(int i=0;i<_q1.n;i++)
    {
        _q1[i]=double(q1[i]);
    }

    c0=c0+b*inv(d,cn/2)*(_q1-c1);
    fpolynome tmp=a-b*inv(d,d.n)*trans(b);
    tmp.n=a.n;
    polynome q0=samplefz(tmp,c0);
    
    return {q0,q1};
}
polynome samplefz(fpolynome f,fpolynome c)
{
    if(f.n<=1) 
        return rand_d(cn/2,0.0,std::min(sqrt(fabs(double(f[0]))),csigma*10));
    else
    {
        fpolynome f0((f.n+1)/2),f1(f.n/2);
        fpolynome _c0((c.n+1)/2),_c1(c.n/2);
//std::cout<<f0<<endl;
        for(int i=0,j=0;i<c.n;i+=2,j++)
        {
            f0[j]=f[i];
            _c0[j]=c[i];
        }
        for(int i=1,j=0;i<c.n;i+=2,j++)
        {
            f1[j]=f[i];
            _c1[j]=c[i];
        }

        std::pair<polynome,polynome> q=sample2z(f0,f1,f0,_c0,_c1);
        polynome ans(q.first.n+q.second.n);
        int j=0;
        for(int i=0;i<q.first.n;i++,j+=2)
        {
            ans[j]=q.first[i];
        }
        j=1;
        for(int i=0;i<q.first.n;i++,j+=2)
        {
            ans[j]=q.first[i];
        }

        return ans;
    }
}

matrix perturb(int n,int q,double s,double alpha,matrix rho,matrix v)
{
    int _m=rho.size()-2;
    double z=1.0/(1.0/(alpha*alpha)-1.0/(s*s));
    fpolynome a(n);
    for(int i=0;i<_m;i++)
    {
        fpolynome tr(n),te(n);
        for(int j=0;j<n;j++)
        {
            tr[j]=double(rho[i][j]);
        }

        a=a+(fpolynome(n)-tr*tr*z);
    }
    a=a+s*s;
    fpolynome b(n);
    for(int i=0;i<_m;i++)
    {
        fpolynome tr(n),te(n);
        for(int j=0;j<n;j++)
        {
            tr[j]=double(rho[i][j]);
        }
        for(int j=0;j<n;j++)
        {
            te[j]=double(v[i][j]);
        }

        a=a+(fpolynome(n)-tr*te*z);
    }
    fpolynome d(n);
    for(int i=0;i<_m;i++)
    {
        fpolynome te(n);
        for(int j=0;j<n;j++)
        {
            te[j]=double(v[i][j]);
        }

        d=d+(fpolynome(n)-te*te*z);
    }
    d=d+s*s;

    matrix _q(_m);
    for(int i=2;i<_m;i++)
    {
        _q.push_back(rand_d(n,0.0,sqrt(s*s-alpha*alpha)));
    }

    polynome tc0(n,q),tc1(n,q);
    for(int i=2;i<_m;i++)
    {
        tc0=tc0+rho[i]*_q[i];
        tc1=tc1+v[i]*_q[i];
    }

    fpolynome c0(n),c1(n);
    double tmp=-(alpha*alpha/(s*s-alpha*alpha));
    for(int i=0;i<n;i++)
    {
        c0[i]=double(tc0[i])*tmp;
        c1[i]=double(tc1[i])*tmp;
    }

    std::pair<polynome,polynome> p=sample2z(a,b,d,c0,c1);
    _q[0]=p.first%q;
    _q[1]=p.second%q;
    for(int i=0;i<_m+2;i++)
    {
        _q[i]=_q[i]%q;
    }

    return _q;
}

matrix sampleg(double s,polynome u,ll q)
{
    int k=ilog2(q)+1;
    double sigma=s/3.0;
    fpolynome l,h,d;
    l[0]=sqrt(2.0*(1.0+1.0/k)+1);
    h[0]=0.0;
    d[0]=double(q%2)/2.0;
    ll tq=q>>1;
    for(int i=1;i<k;i++,tq>>=1)
    {
        l[i]=sqrt(2.0*(1.0+1.0/(k-i)));
        h[i]=sqrt(2.0*(1.0-1.0/(k-i+1)));
        d[i]=(d[i-1]+tq%2)/2;
    }

    matrix Z;
    for(int i=0;i<u.n;i++)
    {
        ll tv=(u[i]+q)%q;

        std::vector<ll> p;
        
        double beta=0.0;

        std::random_device rd;
        std::default_random_engine e(rd());
        std::vector<double> ss,bb,cc;
        std::vector<ll> z;
        for(int i=0;i<l.n;i++)
        {
            cc.push_back(beta/l[i]);
            ss.push_back(sigma/l[i]);
            std::normal_distribution<> nd(cc[i],sigma/l[i]);
            z.push_back(ll(round(nd(e))));
            bb.push_back(-z[i]*h[i]);
        }
        p.push_back(5LL*z[0]+2LL*z[1]);
        for(int i=1;i<k-1;i++)
        {
            p.push_back(2*(z[i-1]+z[i]+z[i+1]));
        }
        p.push_back(2*(z[k-2]+z[k-1]));
        

        fpolynome c(k);
        c[0]=(double(tv%2)-p[0])/2;
        for(int j=1;j<k;j++,tv>>=1)
        {
            c[j]=(c[j-1]+tv%2-p[j])/2;
        }
        
        std::vector<ll> zz(k);
        std::normal_distribution<> nd0(-cc[k-1]/d[k-1],ss[i]/d[k-1]);
        zz[k-1]=ll(round(nd0(e)));
        for(int i=0;i<k-1;i++)
        {
            std::normal_distribution<> nd(sigma,zz[k-1]*d[i]-c[i]);
            zz[i]=ll(round(nd(e)));
        }

        polynome t(k);
        tv=(u[i]+q)%q;
        std::vector<ll> ttv;
        for(int i=0;i<k;i++)
        {
            ttv.push_back(tv%2);
            tv>>=1;
        }
        t[0]=2*z[0]+q%2*z[k-1]+ttv[0];
        q>>=1;
        for(int j=1;j<k-1;j++,q>>=1)
        {
            t[j]=2*z[j]-z[j-1]+q%2*z[k-1]+ttv[j];
        }
        t[k-1]=q%2*z[k-1]-z[k-2]+ttv[k-1];
        Z.push_back(t);
    }

    matrix ans;
    for(int i=0;i<k;i++)
    {
        polynome tans(u.n);
        for(int j=0;j<u.n;j++)
        {
            tans[j]=Z[j][i];
        }
        ans.push_back(tans%q);
    }

    return ans;
}

matrix gausssamp(matrix a,matrix rho,matrix v,polynome beta,double sigma)
{
    // matrix p=perturb(cn/2,cq,1.3*3.0*sigma*sigma*(sqrt(cn*rho.size()/2)+sqrt(cn)+4.7),3.0*sigma,rho,v);
    // for(int i=0;i<a.size();i++)
    // {
    //     beta=beta-a[i]*p[i];
    // }
    matrix z=sampleg(sigma,beta,cq);

    polynome ez(cn/2,cq),rz(cn/2,cq);
    for(int i=0;i<v.size();i++)
    {
        ez=ez+v[i]*z[i];
        rz=rz+rho[i]*z[i];
    }

    matrix x;
    x.push_back(ez);
    x.push_back(rz);
    for(int i=2;i<rho.size();i++)
    {
        x.push_back(z[i-2]);
    }

    return x;
}

#endif
