// Exact quotient of balanced hexaploid phase histories by homolog relabelling.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;
using Clock=chrono::steady_clock;
uint64_t encode(array<unsigned char,6> a){sort(a.begin(),a.end());uint64_t k=0;for(int h=0;h<6;++h)k|=uint64_t(a[h])<<(8*h);return k;}
array<unsigned char,6> decode(uint64_t k){array<unsigned char,6>a{};for(int h=0;h<6;++h)a[h]=(k>>(8*h))&255;return a;}
int main(int argc,char**argv){
 int d=argc>1?atoi(argv[1]):6;string file=argc>2?argv[2]:"graph.bin";
 if(d<1||d>7)return 2;
 auto start=Clock::now();vector<int>masks;for(int m=0;m<64;++m)if(__builtin_popcount(unsigned(m))==3)masks.push_back(m);
 int cap=(1<<d)-1;array<unsigned char,6>a{0,0,0,(unsigned char)cap,(unsigned char)cap,(unsigned char)cap};
 vector<uint64_t>keys{encode(a)};vector<uint32_t>edges;unordered_map<uint64_t,uint32_t>index;
 if(d==7){keys.reserve(2000000);edges.reserve(40000000);index.reserve(2000000);}else{keys.reserve(120000);edges.reserve(2400000);index.reserve(120000);}
 index[keys[0]]=0;
 for(size_t s=0;s<keys.size();++s){
  auto old=decode(keys[s]);
  for(int m:masks){array<unsigned char,6>b{};for(int h=0;h<6;++h)b[h]=((old[h]<<1)|((m>>h)&1))&cap;
   uint64_t k=encode(b);auto got=index.find(k);uint32_t id;
   if(got==index.end()){id=keys.size();index[k]=id;keys.push_back(k);}else id=got->second;
   edges.push_back(id);
  }
  if(s&&s%200000==0)cerr<<"visited "<<s<<" discovered "<<keys.size()<<" seconds "<<chrono::duration<double>(Clock::now()-start).count()<<"\n";
 }
 uint64_t labelled=0;vector<double>prior;
 for(auto key:keys){auto r=decode(key);int denominator=1;for(int h=0;h<6;){int j=h+1;while(j<6&&r[j]==r[h])++j;for(int z=2;z<=j-h;++z)denominator*=z;h=j;}int multiplicity=720/denominator;labelled+=multiplicity;prior.push_back(multiplicity);}
 uint64_t expected=1;for(int j=0;j<d;++j)expected*=20;
 if(labelled!=expected){cerr<<"Incomplete orbit space "<<labelled<<" vs "<<expected<<"\n";return 3;}
 for(auto&x:prior)x/=double(expected);
 // Every labelled column is uniform; orbit weights must remain stationary.
 vector<double>next(keys.size(),0);for(size_t s=0;s<keys.size();++s)for(int m=0;m<20;++m)next[edges[20*s+m]]+=prior[s]/20;
 double err=0;for(size_t s=0;s<keys.size();++s)err=max(err,abs(next[s]-prior[s]));
 ofstream out(file,ios::binary);uint64_t n=keys.size();out.write((char*)&d,sizeof(d));out.write((char*)&n,sizeof(n));out.write((char*)keys.data(),n*sizeof(uint64_t));out.write((char*)edges.data(),n*20*sizeof(uint32_t));out.write((char*)prior.data(),n*sizeof(double));
 cout<<"{\"history\":"<<d<<",\"orbits\":"<<n<<",\"labelled_histories\":"<<labelled<<",\"stationarity_error\":"<<err<<",\"seconds\":"<<chrono::duration<double>(Clock::now()-start).count()<<"}\n";
}
