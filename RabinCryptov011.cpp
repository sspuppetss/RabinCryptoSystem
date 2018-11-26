#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <cstring>
#include "mini-gmp.h"
#include "mini-gmp.c"


using namespace std;

struct twolong
{
    long x;
    long y;
};

struct fourlong
{
    long w;
    long x;
    long y;
    long z;
};

struct twovec
{
    vector<int> a;
    vector<int> b;
};

struct fourvec
{
    vector<int> a;
    vector<int> b;
    vector<int> c;
    vector<int> d;
};

long pangkat(long x,long y)
{
    long res;
    res=1;
    for(int i = 0; i<y;i++){
        res = res*x;
    }
    return res;
}

void printvec(vector<int> v)
{
    copy(v.begin(), v.end(), ostream_iterator<int>(cout, ""));
}

vector<int> dectobin(long m)
{
    vector<int>l;
    long temp;
        while(m > 0)
         {
          	temp = m%2;
        	l.insert(l.begin(),temp);
        	m = m/2;
         }
    return l;
}

long bintodec(vector<int> m)
{
    long power,decimal=0;
        for(int i = m.size()-1;i>=0;i--)
         {
		    power = m.size()-i-1;
		    decimal += m[i]*pangkat(2,power);
         }
    return decimal;
}

vector<int> doublebin(vector<int> m)
{
    vector<int>l;
    for(int i = 0; i < m.size(); i++)
    {
       l.push_back(m[i]);
    }
    for(int i = 0; i < m.size(); i++)
    {
       l.push_back(m[i]);
    }
    return l;
}

long kongruen(long m, long n)
{
    long k;
    k = (m-(m%n))/n;
    return k;
}

long enkripsi(long m, long n)
{
    long c;
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_ui(tmp,m);
    mpz_pow_ui(tmp,tmp,2);
    mpz_mod_ui(tmp,tmp,n);
    c = mpz_get_ui(tmp);
    mpz_clear(tmp);
    return c;
}

twolong euclidex(long p, long q)
{
            long fx,fy,tempp,tempq;
            tempp = p;
            tempq = q;
	        long x = 0, y = 1;
	        long lastx = 1, lasty = 0, temp;
	        while (tempq != 0)
	        {
	            long z1 = tempp / tempq;
	            long z2 = tempp % tempq;
	            tempp = tempq;
	            tempq = z2;

	            temp = x;
	            x = lastx - z1 * x;
	            lastx = temp;

	            temp = y;
	            y = lasty - z1 * y;
	            lasty = temp;
	        }
	        fx = lastx;
	        fy = lasty;
            twolong euxy = {fx,fy};
            return euxy;
}

twolong quampq(long c, long p, long q)
{
    twolong mpq,pq;
    pq.x = (p+1)/4;
    pq.y = (q+1)/4;
    mpz_t mp,mp1,mp2;
    mpz_init(mp);
    mpz_init(mp1);
    mpz_init(mp2);
    mpz_set_ui(mp,c);
    mpz_pow_ui(mp1,mp,pq.x);
    mpz_pow_ui(mp2,mp,pq.y);
    mpz_mod_ui(mp1,mp1,p);
    mpz_mod_ui(mp2,mp2,q);
    mpq.x = mpz_get_ui(mp1);
    mpq.y = mpz_get_ui(mp2);
    mpz_clear(mp);
    mpz_clear(mp1);
    mpz_clear(mp2);
    return mpq;
}

fourlong crt(twolong eu, twolong mpq, long p, long q,long n)
{
    fourlong cr;
    cr.w = (eu.x*p* mpq.y + eu.y * q* mpq.x) % n ;
    cr.x = (-cr.w)%n;
    cr.y = (eu.y*q* mpq.x - eu.x* p* mpq.y ) % n ;
    cr.z = (-cr.y)%n;
    if(cr.w<0)cr.w=n+cr.w;
    if(cr.x<0)cr.x=n+cr.x;
    if(cr.y<0)cr.y=n+cr.y;
    if(cr.z<0)cr.z=n+cr.z;
    return cr;
}

fourlong krstu(long k, fourlong rstu,long n)
{
    fourlong l;
    l.w = (k*n)+rstu.w;
    l.x = (k*n)+rstu.x;
    l.y = (k*n)+rstu.y;
    l.z = (k*n)+rstu.z;
    return l;
}

long dekrip( vector <int> c)
{
    vector <int> l1,l2;
    long res = 0,temp;
    if(c.size()%2==0)
    {
        temp = c.size()/2;
    }else{
        temp = (c.size()+1)/2;
    }
    for(int i=0;i<temp;i++){
        l1.push_back(c[i]);
    }
    for(int j=c.size()-1;j>=temp;j--){
        l2.insert(l2.begin(),c[j]);
    }
    int k;
    for (k = 0; k<temp;k++){
        if(l1[k]!=l2[k]){
            break;
        }
    }
    if(k==temp){
        res = bintodec(l1);
    }
    return res;
}

int main()
{
    long m,p,q,n,mbaru,k,c;
    vector<int> mbin,mbin2;
    twolong eu,mpq,test;
    fourlong rstu,nrstu,cm;
    fourvec cbin;

    p = 23;
    q = 11;
    n = p * q;
    m = 8;

    //enkripsi
    mbin = dectobin(m);
    mbin2 = doublebin(mbin);
    mbaru = bintodec(mbin2);
    k = kongruen(mbaru, n);
    c = enkripsi(mbaru,n);

    cout<<"=============== Proses Enkripsi ============="<<endl;
    cout<<"Pesan / Plain Text : "<<m<<endl;
    cout<<endl<<"Kunci Publik : "<<n<<endl;
    cout<<"Biner Pesan: ";
    printvec(mbin);cout<<endl;
    cout<<"Duplikasi Biner pesan: ";
    printvec(mbin2);cout<<endl;
    cout<<"Pesan Baru : "<<mbaru<<endl;
    cout<<"Kongruen : "<<k<<endl;
    cout<<endl<<"Chiper Text: "<<c<<endl;

    //decripsi
    eu = euclidex(p,q);
    mpq = quampq(c,p,q);
    rstu = crt(eu,mpq,p,q,n);
    nrstu = krstu(k,rstu,n);
    cbin.a = dectobin(nrstu.w);
    cbin.b = dectobin(nrstu.x);
    cbin.c = dectobin(nrstu.y);
    cbin.d = dectobin(nrstu.z);

    cout<<endl<<"=============== Proses Dekripsi ============="<<endl;
    cout<<"============= EuclidEx Algoritma ============"<<endl;
    cout<<"Yp : "<<eu.x<<endl;
    cout<<"Yq : "<<eu.y<<endl;
    cout<<endl<<"=== Nilai Chiper Kuadrat Terhadap p dan q ==="<<endl;
    cout<<"mp : "<<mpq.x<<endl;
    cout<<"mq : "<<mpq.y<<endl;
    cout<<endl<<"=== Chinese Remainder Theorem + Kongruen ===="<<endl;
    cout<<"=============== Hasil Dekripsi =============="<<endl;
    cout<<"r : "<<nrstu.w<<"\t ==> ";
    printvec(cbin.a);cout<<endl;
    cout<<"s : "<<nrstu.x<<"\t ==> ";
    printvec(cbin.b);cout<<endl;
    cout<<"t : "<<nrstu.y<<"\t ==> ";
    printvec(cbin.c);cout<<endl;
    cout<<"u : "<<nrstu.z<<"\t ==> ";
    printvec(cbin.d);cout<<endl;

    cout<<endl<<"========================="<<endl;
    cm.w = dekrip(cbin.a);
    cm.x = dekrip(cbin.b);
    cm.y = dekrip(cbin.c);
    cm.z = dekrip(cbin.d);
    if(cm.w!=0){
        cout<<"Pesan Asli : "<<cm.w<<endl;
    }else if(cm.x!=0){
        cout<<"Pesan Asli : "<<cm.x<<endl;
    }else if(cm.y!=0){
        cout<<"Pesan Asli : "<<cm.y<<endl;
    }else if(cm.z!=0){
        cout<<"Pesan Asli : "<<cm.z<<endl;
    }else {
        cout<<"Error"<<endl;
    }
    cout<<"========================="<<endl;
}
