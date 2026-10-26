// Exact finite-state marginalization. Checkpoints bound resident memory.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=5){cerr<<"graph tables checkpoints output\n";return 2;}
 auto start=chrono::steady_clock::now();ifstream graph(argv[1],ios::binary);int d;uint64_t n;
 graph.read((char*)&d,4);graph.read((char*)&n,8);if(!graph||d<1||d>7)return 3;
 vector<uint64_t>keys(n);vector<uint32_t>edges(n*20);vector<double>alpha(n),next(n),beta(n,1),previous(n);
 graph.read((char*)keys.data(),n*8);graph.read((char*)edges.data(),n*20*4);graph.read((char*)alpha.data(),n*8);if(!graph)return 4;
 ifstream input(argv[2],ios::binary);int stages,reads,patterns;input.read((char*)&stages,4);input.read((char*)&reads,4);input.read((char*)&patterns,4);
 if(patterns!=(1<<(d+1))||reads!=2)return 5;
 vector<double>table(size_t(stages)*reads*3*patterns);input.read((char*)table.data(),table.size()*8);if(!input)return 6;
 vector<int>masks;for(int m=0;m<64;++m)if(__builtin_popcount((unsigned)m)==3)masks.push_back(m);
 vector<double>scales(stages),answers(size_t(stages)*reads*2);
 fstream checkpoint(argv[3],ios::binary|ios::in|ios::out|ios::trunc);double logZ=0;
 auto sums=[&](uint64_t key,int mask,const double*t,double*tot){
  for(int j=0;j<6;++j)tot[j]=0;
  for(int h=0;h<6;++h){int pat=2*((key>>(8*h))&255)+((mask>>h)&1);
   for(int r=0;r<2;++r)for(int a=0;a<3;++a)tot[3*r+a]+=t[(3*r+a)*patterns+pat];}
 };
 for(int t=0;t<stages;++t){
  checkpoint.write((char*)alpha.data(),n*8);fill(next.begin(),next.end(),0);const double*tab=table.data()+size_t(t)*6*patterns;
  for(uint64_t s=0;s<n;++s)if(alpha[s]){
   double z[6];for(int j=0;j<20;++j){sums(keys[s],masks[j],tab,z);double L=z[0]*z[3]/36;
    next[edges[20*s+j]]+=alpha[s]*L/20;}
  }
  long double norm=0;for(auto v:next)norm+=v;scales[t]=double(norm);if(!(norm>0))return 7;
  logZ+=log(double(norm));for(auto&v:next)v/=double(norm);alpha.swap(next);
  if(t%8==0)cerr<<"forward "<<t<<" logZ "<<setprecision(12)<<logZ<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
 }
 checkpoint.flush();
 for(int t=stages-1;t>=0;--t){
  checkpoint.seekg(uint64_t(t)*n*8);checkpoint.read((char*)alpha.data(),n*8);if(!checkpoint)return 8;
  const double*tab=table.data()+size_t(t)*6*patterns;long double mass=0,p0[2]={0,0},errors[2]={0,0};
  for(uint64_t s=0;s<n;++s){double z[6],b=0;
   for(int j=0;j<20;++j){sums(keys[s],masks[j],tab,z);double L=z[0]*z[3]/36;
    double x=L*beta[edges[20*s+j]]/(20*scales[t]);b+=x;double w=alpha[s]*x;mass+=w;
    for(int r=0;r<2;++r){p0[r]+=w*z[3*r+1]/z[3*r];errors[r]+=w*z[3*r+2]/z[3*r];}
   }previous[s]=b;
  }
  if(abs(double(mass)-1)>1e-8){cerr<<"posterior mass "<<double(mass)<<"\n";return 9;}
  for(int r=0;r<2;++r){answers[(t*2+r)*2]=double(p0[r]);answers[(t*2+r)*2+1]=double(errors[r]);}
  beta.swap(previous);
  if(t%8==0)cerr<<"backward "<<t<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
 }
 ofstream output(argv[4],ios::binary);output.write((char*)&logZ,8);output.write((char*)answers.data(),answers.size()*8);
 cout<<setprecision(16)<<"{\"log_evidence\":"<<logZ<<",\"reads\":"<<stages*2<<",\"states\":"<<n<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
