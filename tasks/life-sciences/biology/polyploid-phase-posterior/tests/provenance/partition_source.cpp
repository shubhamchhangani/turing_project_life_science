// Independent author expectation source: ratios of separately evaluated
// partition functions, with no backward messages or reference posterior output.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=4)return 2;
 ifstream g(argv[1],ios::binary);int d;uint64_t N;g.read((char*)&d,4);g.read((char*)&N,8);
 vector<uint64_t>key(N);vector<uint32_t>dest(20*N);vector<double>prior(N);
 g.read((char*)key.data(),8*N);g.read((char*)dest.data(),80*N);g.read((char*)prior.data(),8*N);if(!g)return 3;
 ifstream f(argv[2],ios::binary);int T,R,P;f.read((char*)&T,4);f.read((char*)&R,4);f.read((char*)&P,4);
 vector<double>lookup(size_t(T)*6*P);f.read((char*)lookup.data(),lookup.size()*8);if(!f||R!=2)return 4;
 vector<int>choice;for(int x=0;x<64;++x)if(__builtin_popcount((unsigned)x)==3)choice.push_back(x);
 vector<double>x=prior,y(N),prefix=prior;vector<double>prefix_log(T+1);double baseline;
 auto advance=[&](vector<double>&v,vector<double>&out,int t,int changed_read,int statistic){
  fill(out.begin(),out.end(),0.0);const double*L=lookup.data()+size_t(t)*6*P;
  int component0=(changed_read==0?statistic:0),component1=(changed_read==1?statistic:0);
  const double*A=L+component0*P;const double*B=L+(3+component1)*P;
  for(uint64_t i=0;i<N;++i){
   unsigned code[6];for(int h=0;h<6;++h)code[h]=2*((key[i]>>(8*h))&255);
   double weight=v[i]/720.0; // 20 genotype choices and two independent 1/6 origins.
   for(int j=0;j<20;++j){double a=0,b=0;unsigned c=choice[j];
    for(int h=0;h<6;++h){unsigned p=code[h]+((c>>h)&1);a+=A[p];b+=B[p];}
    out[dest[20*i+j]]+=weight*a*b;
   }
  }
  long double z=0;for(double a:out)z+=a;
  if(!(z>0))return -HUGE_VAL;
  double scale=double(z);for(double&a:out)a/=scale;v.swap(out);return log(scale);
 };
 auto begin=chrono::steady_clock::now();
 for(int t=0;t<T;++t)prefix_log[t+1]=prefix_log[t]+advance(x,y,t,-1,0);
 baseline=prefix_log.back();ofstream out(argv[3]);out<<setprecision(17)<<"{\"log_evidence\":"<<baseline<<",\"read_expectations\":[";
 bool first=true;
 for(int t=0;t<T;++t){
  for(int r=0;r<2;++r){double result[2];
   for(int metric=1;metric<=2;++metric){x=prefix;double value=prefix_log[t];
    for(int j=t;j<T;++j)value+=advance(x,y,j,j==t?r:-1,j==t?metric:0);
    result[metric-1]=exp(value-baseline);
   }
   if(!first)out<<",";first=false;out<<"["<<result[0]<<","<<result[1]<<"]";out.flush();
  }
  advance(prefix,y,t,-1,0);
  cerr<<"certified stage "<<t<<" / "<<T<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-begin).count()<<"\n";
 }
 out<<"],\"method\":\"independent modified-partition ratios\",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-begin).count()<<"}\n";
}
