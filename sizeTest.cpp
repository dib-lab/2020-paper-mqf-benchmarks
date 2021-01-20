#include "gqf.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>
#include<iostream>
#include <fstream>
#include "hashutil.h"
#include "kmer.h"
#include <vector>
#include <chrono>
#include<cmath>
#include <random>
#include <algorithm>
#include "Benchmark_Utils_Script/generators.hpp"
#include "Benchmark_Utils_Script/countingStructure.hpp"
#include "CLI11.hpp"
#include "countmin-khmer/storage.hh"

using namespace std;

#include <iostream>

inline uint64_t estimateMemory(uint64_t nslots, uint64_t slotSize, uint64_t fcounter, uint64_t tagSize) {
    uint64_t SLOTS_PER_BLOCK_t = 64;
    uint64_t xnslots = nslots + 10 * sqrt((double) nslots);
    uint64_t nblocks = (xnslots + SLOTS_PER_BLOCK_t - 1) / SLOTS_PER_BLOCK_t;
    uint64_t blocksize = 17;

    return ((nblocks) * (blocksize + 8 * (slotSize + fcounter + tagSize))) / 1024;

}

bool isEnough(vector<uint64_t>& histogram, uint64_t noSlots, uint64_t fixedSizeCounter, uint64_t slotSize) {
    // cout<<"noSlots= "<<noSlots<<endl
    //     <<"fcounter= "<<fixedSizeCounter<<endl
    //     <<"slot size= "<<numHashBits<<endl;

    noSlots = (uint64_t) ((double) noSlots * 0.95);
    for (uint64_t i = 1; i < 1000; i++) {
        uint64_t usedSlots = 1;

        if (i > ((1ULL) << fixedSizeCounter) - 1) {
            uint64_t nSlots2 = 0;
            __uint128_t capacity;
            do {
                nSlots2++;
                capacity = ((__uint128_t) (1ULL) << (nSlots2 * slotSize + fixedSizeCounter)) - 1;
                //  cout<<"slots num "<<nSlots2<<" "<<capacity<<endl;
            } while ((__uint128_t) i > capacity);
            usedSlots += nSlots2;
        }
        //cout<<"i= "<<i<<"->"<<usedSlots<<" * "<<histogram[i]<<endl;
        if (noSlots >= (usedSlots * histogram[i])) {
            noSlots -= (usedSlots * histogram[i]);
        } else {
            //  cout<<"failed"<<endl<<endl;
            return false;
        }

    }
    //cout<<"success"<<endl<<endl;
    return true;
}
void estimateMemRequirement(vector<uint64_t>& histogram,
                            uint64_t numHashBits, uint64_t tagSize,
                            uint64_t *res_noSlots, uint64_t *res_fixedSizeCounter, uint64_t *res_memory) {
    uint64_t noDistinctKmers = 0, totalNumKmers=0;
    *res_memory = numeric_limits<uint64_t>::max();
    for (int i = 8; i < 64; i++) {
        uint64_t noSlots = (1ULL) << i;
        if (noSlots < noDistinctKmers)
            continue;
        bool moreWork = false;
        uint64_t slotSize = numHashBits - log2((double) noSlots);
        for (uint64_t fixedSizeCounter = 1; fixedSizeCounter < slotSize; fixedSizeCounter++) {
            if (isEnough(histogram, noSlots, fixedSizeCounter, slotSize)) {
                uint64_t tmpMem = estimateMemory(noSlots, slotSize, fixedSizeCounter, tagSize);
                if (*res_memory > tmpMem) {
                    *res_memory = tmpMem;
                    *res_fixedSizeCounter = fixedSizeCounter;
                    *res_noSlots = i;
                    moreWork = true;
                } else {
                    break;
                }
            }

        }
        if (!moreWork && *res_memory != numeric_limits<uint64_t>::max())
            break;
    }
    if (*res_memory == numeric_limits<uint64_t>::max()) {
        throw std::overflow_error(
                "Data limits exceeds MQF capabilities(> uint64). Check if ntcard file is corrupted");
    }


}

int main(int argc, char const *argv[]) {

  cout<<"Command = ";
  for(int i=0;i<argc;i++)
  {
    cout<<argv[i]<<" ";
  }
  cout<<endl;

  CLI::App app;
  string distribution="zipfian";
  double zipifian_coefficient=3.0;
  double uniform_coefficient=5;
  string dataStrucureInput="mqf";
  double fpr=0.01;
  uint64_t countErrorMargin=5;
  uint64_t kSize=25;
  string fqPath="";

  uint64_t qbits=25;


  uint64_t num_elements=10000000;


  app.add_option("-d,--distribution",distribution,
  "Distributions which items will be drew from. options are zipfian,kmers. Default is zipfian");

  app.add_option("-z",zipifian_coefficient,"Zipifian Coeffiecient. Default =3")->group("Zipfian Distribution Options");
  app.add_option("-u",uniform_coefficient,
  "Frequency of items in the uinform distribution. Default =5")
  ->group("Uniform Distribution Options");
  app.add_option("-k",fqPath,"Fastq file to draw the kmers");

  // app.add_option("-s,--datastructure",dataStrucureInput,
  // "Datastructure to be benchmarked. Options are mqf, lmqf, cqf, and countmin. Default is mqf");

  // app.add_option("-s,--size",qbits
  // ,"No slots equals to 2^(s). default is 25")
  // ->group("Counting Strucutre Options");

  app.add_option("-f,--fpr",fpr
  ,"False positive rate. default is 0.01")
  ->group("Counting Strucutre Options");

  app.add_option("-n,--num_elements",num_elements,
  "Number of elements to generate. Default is 10M");

  app.add_option("-c,--countsError",countErrorMargin
  ,"Counts error margin. default is 5")
  ->group("Countmin options");





  CLI11_PARSE(app, argc, (char**)argv);

  uint64_t p=(uint64_t)log2((double)num_elements/fpr);

  cout<<"Num Hashbits= "<<p<<endl;
  uint64_t estimatedQ=log2((double)num_elements*1.1)+1;// worst case zipifian distribution
  cout<<"Estimated Q= "<<estimatedQ<<endl;

  vector<countingStructure*> dataStructures;

  dataStructures.push_back(new MQF(estimatedQ+1,p-estimatedQ-1,1));
  dataStructures.push_back(new MQF(estimatedQ,p-estimatedQ,5));
  dataStructures.push_back(new MQF(estimatedQ,p-estimatedQ,4));
  dataStructures.push_back(new MQF(estimatedQ,p-estimatedQ,3));
  dataStructures.push_back(new MQF(estimatedQ,p-estimatedQ,2));
  dataStructures.push_back(new MQF(estimatedQ,p-estimatedQ,1));
  dataStructures.push_back(new MQF(estimatedQ-1,p-estimatedQ+1,5));
  dataStructures.push_back(new MQF(estimatedQ-1,p-estimatedQ+1,4));
  dataStructures.push_back(new MQF(estimatedQ-1,p-estimatedQ+1,3));
  dataStructures.push_back(new MQF(estimatedQ-1,p-estimatedQ+1,2));
  dataStructures.push_back(new MQF(estimatedQ-1,p-estimatedQ+1,1));
  dataStructures.push_back(new MQF(estimatedQ-2,p-estimatedQ+2,5));
  dataStructures.push_back(new MQF(estimatedQ-2,p-estimatedQ+2,4));
  dataStructures.push_back(new MQF(estimatedQ-2,p-estimatedQ+2,3));
  dataStructures.push_back(new MQF(estimatedQ-2,p-estimatedQ+2,2));
  dataStructures.push_back(new MQF(estimatedQ-2,p-estimatedQ+2,1));


  dataStructures.push_back(new CQF(estimatedQ+1,p-estimatedQ-1));
  dataStructures.push_back(new CQF(estimatedQ,p-estimatedQ));
  dataStructures.push_back(new CQF(estimatedQ-1,p-estimatedQ+1));
   dataStructures.push_back(new CQF(estimatedQ-2,p-estimatedQ+2));

  //dataStructures.push_back(new LMQF(estimatedQ-1,estimatedQ-2,p-estimatedQ+1,2));
  //dataStructures.push_back(new LMQF(estimatedQ-1,estimatedQ-3,p-estimatedQ+1,2));

  uint64_t BufferSize=num_elements/25;
  uint64_t num_queries=100;


  int countMinDepth=log(1.0/(fpr))+1;
  double e=2.71828;
//  uint64_t countMinWidth =e*(double)num_elements/(double) (countErrorMargin+1) ;
  uint64_t countMinWidth=dataStructures[0]->size/countMinDepth;

  //cout<<"Count Min Sketch Width= "<<countMinWidth<<endl;
  //cout<<"Count Min Sketch Depth= "<<countMinDepth<<endl;

  //dataStructures.push_back(new countminKhmer((uint64_t)((double)countMinWidth*1.2),countMinDepth));
//  dataStructures.push_back(new countminKhmer((uint64_t)((double)countMinWidth*1),countMinDepth));
//  dataStructures.push_back(new countminKhmer((uint64_t)((double)countMinWidth*0.8),countMinDepth));


  //dataStructures.push_back(new countmin(countMinWidth,countMinDepth));


  uint64_t range=(1ULL<<(int)p);
  srand (1);

  map<uint64_t,uint64_t> gold;


 
  generator *g;
  if(distribution=="uniform")
  {
    g=new uniformGenerator(num_elements,range,uniform_coefficient);
  }
  else if(distribution=="kmers")
  {
    g= new kmersGenerator(num_elements,fqPath,kSize);
  }
  else if(distribution=="zipfian"){
    g= new zipfGenerator(num_elements,range,zipifian_coefficient);
  }
  else{
    cout<<"Unknown Distribution"<<endl;
    return -1;
  }




  //
  //
  uint64_t countedKmers=0;
  vector<uint64_t> insertions(BufferSize);
  //set<uint64_t> succesfullQueries(num_queries);
  uint64_t queryTop=0;
  // auto now = std::chrono::high_resolution_clock::now();
  // auto prev=now;
  // now = std::chrono::high_resolution_clock::now();
  //
  //
  //
  // vector<pair<int,int> > lmqfSpace;
  // vector<pair<double,double> > result;
  // vector<pair<double,double> > results_slotsUsed;
  // cout<<"Capacity\t1M insertions per second"<<endl;
  // double capacity=dataStructure->space();
  //
  //
  // uint64_t old_capacity=capacity;
  // uint64_t insertiontime=0;
  // uint64_t querytime=0;
  // uint64_t curr_item;
  auto now=std::chrono::high_resolution_clock::now();
  auto prev=std::chrono::high_resolution_clock::now();
  auto microseconds = (chrono::duration_cast<chrono::microseconds>(now-prev)).count();
  uint64_t curr_item;

// check the size of cqf
   while(g->hasMore())
   {
     //     cerr<<"Number of unique items generated so far ="<<g->nunique_items<<endl;
     uint64_t tmp_inserted=0;
     for(int j=0;j<BufferSize;j++){
       if(!g->getElement(curr_item))
         break;
       tmp_inserted++;
    //   cout<<"sss"<<curr_item<<"\n";
       insertions[j]=curr_item;
    //   gold[curr_item]++;
       // if(queryTop<num_queries)
       //    succesfullQueries[queryTop++]=curr_item;
     }
     //    cerr<<"Finished generation"<<endl;
     countedKmers+=BufferSize;
     for(auto structure: dataStructures){
       if(structure->space()>=95)
       {
         structure->failed=true;
         continue;
       }
       //  cerr<<"Inserting "<<BufferSize<<" to "<<structure->name<<"("<<structure->space()<<")"<<endl;
       prev=std::chrono::high_resolution_clock::now();

       for(int j=0;j<tmp_inserted;j++){
          structure->insert(insertions[j],1);
          if(structure->space()>=95)
          {
            structure->failed=true;
            break;
          }

        }
       now=std::chrono::high_resolution_clock::now();
       microseconds = (chrono::duration_cast<chrono::microseconds>(now-prev)).count();
       microseconds/=1000;
       structure->insertionTime+=microseconds;
     }
  }
  cout<<"Final Number of unique items generated="<<g->nunique_items<<endl;

  // for(auto structure: dataStructures){
  //   if(structure->failed)
  //   {
  //     structure->fpr=1;
  //     structure->fpr5=1;
  //     structure->fpr10=1;
  //     structure->fpr20=1;
  //     structure->fpr30=1;
  //     continue;
  //
  //   }
  //   cerr<<"Querying  from "<<structure->name<<endl;
  //   prev=std::chrono::high_resolution_clock::now();
  //
  //   for(auto a :g->newItems)
  //   {
  //      uint64_t tmpCount=structure->query(a);
  //
  //      if(tmpCount>0)
  //        structure->fpr++;
  //      if(tmpCount>5)
  //        structure->fpr5++;
  //      if(tmpCount>10)
  //        structure->fpr10++;
  //      if(tmpCount>20)
  //         structure->fpr20++;
  //      if(tmpCount>30)
  //         structure->fpr30++;
  //   }
  //   uint64_t numQueries5=0,
  //            numQueries10=0,
  //            numQueries20=0,
  //            numQueries30=0;
  //
  //   // for(auto it=gold.begin();it!=gold.end();it++)
  //   // {
  //   //   uint64_t tmpCount=structure->query(it->first);
  //   //   if(it->second<5){
  //   //     numQueries5++;
  //   //     if(tmpCount>5)
  //   //       structure->fpr5++;
  //   //   }
  //   //   if(it->second<10){
  //   //     numQueries10++;
  //   //     if(tmpCount>10)
  //   //       structure->fpr10++;
  //   //   }
  //   //   if(it->second<20){
  //   //     numQueries20++;
  //   //     if(tmpCount>20)
  //   //       structure->fpr20++;
  //   //   }
  //   //   if(it->second<30){
  //   //     numQueries30++;
  //   //     if(tmpCount>30)
  //   //       structure->fpr30++;
  //   //   }
  //   //
  //   // }
  //   // now=std::chrono::high_resolution_clock::now();
  //   // structure->fpr/=double(g->newItems.size());
  //   // structure->fpr5/=double(g->newItems.size()+numQueries5);
  //   // structure->fpr10/=double(g->newItems.size()+numQueries10);
  //   // structure->fpr20/=double(g->newItems.size()+numQueries20);
  //   // structure->fpr30/=double(g->newItems.size()+numQueries30);
  //
  //   microseconds = (chrono::duration_cast<chrono::microseconds>(now-prev)).count();
  //   microseconds/=1000;
  //   structure->queryTime+=microseconds;
  //
  //
  // }


  cout<<"Number of insertions = "<<countedKmers<<endl;
//  cout<<"Number of succesfull lookups = "<<num_queries<<endl;
//  cout<<"Number of non succesfull lookups = "<<num_queries<<endl;


  cout<<"Name"<<"\t"
      <<"Parametes"<<"\t"
      <<"Size"<<"\t"
      <<"Space"<<"\t"
      <<"Succeded"<<"\t"
      <<"fpr"<<"\t"
      <<"fpr5"<<"\t"
      <<"fpr10"<<"\t"
      <<"fpr20"<<"\t"
      <<"fpr30"<<endl;

  for(auto structure: dataStructures){
      cout<<structure->name<<"\t"
          <<structure->parameters<<"\t"
          <<structure->size/(1024*1024)<<"MB\t"
          <<structure->space()<<"\t"
          <<!structure->failed<<"\t"
          <<structure->fpr<<"\t"
          <<structure->fpr5<<"\t"
          <<structure->fpr10<<"\t"
          <<structure->fpr20<<"\t"
          <<structure->fpr30<<endl;
  }
  vector<uint64_t> histogram(10000);
  for(auto structure: dataStructures)
  {
      if(structure->name=="MQF" && !structure->failed)
      {
          QF* mqf=((MQF*)structure)->get_MQF();
          QFi iterator;
          qf_iterator(mqf,&iterator,0);
          while(!qfi_end(&iterator))
          {
              uint64_t key,count,tag;
              qfi_get(&iterator,&key,&tag,&count);
              histogram[count]++;
              qfi_next(&iterator);
          }
          uint64_t res_noSlots,res_fixedSizeCounter, res_memory;
          estimateMemRequirement(histogram,p,0,&res_noSlots, &res_fixedSizeCounter,&res_memory);
          cout<<"Estimated Q = "<<res_noSlots<<endl;
          cout<<"Estimated Fixed Size counter = "<<res_fixedSizeCounter<<endl;
          cout<<"Estimated memory = "<<res_memory<<endl;

          break;
      }

  }
  //
  //   now = std::chrono::high_resolution_clock::now();
  //   if(dataStrucureInput=="bmqf")
  //   {
  //
  //   }
  //   else{
  //     for(int j=0;j<BufferSize;j++)
  //     {
  //       dataStructure->query(input[j]);
  //     }
  //   }
  //
  //   prev=now;
  //   now = std::chrono::high_resolution_clock::now();
  //
  //   microseconds = (chrono::duration_cast<chrono::microseconds>(now-prev)).count();
  //   querytime=microseconds;
  //   {
  //     double millionInsertionSecond=(double)countedKmers/(double)insertiontime;
  //     cout<<capacity<<"\t"<<g->get_generated_elements()<<"\t"<<millionInsertionSecond;
  //     if (dataStrucureInput=="lmqf"){
  //       pair<int, int> tmp;
  //       layeredMQF* lmqf=((LMQF*)dataStructure)->get_MQF();
  //       tmp.first=qf_space(lmqf->firstLayer_singletons);
  //       tmp.second=qf_space(lmqf->secondLayer);
  //       int64_t tmp2=(int64_t)slotsUsedInCounting(lmqf->firstLayer_singletons);
  //       tmp2-=(int64_t)lmqf->firstLayer_singletons->metadata->noccupied_slots;
  //       cout<<"\t"<<tmp.first<<"\t"<<tmp.second<<"\t"<<tmp2;
  //     }
  //     cout<<endl;
  //     old_capacity=capacity;
  //
  //     capacity=dataStructure->space();
  //
  //     double millionQuerysSecond=(double)countedKmers/(double)querytime;
  //     result.push_back(make_pair(capacity,(double)millionQuerysSecond));
  //   }
  //
  //   capacity=dataStructure->space();
  //   if(dataStrucureInput=="mqf")
  //   {
  //     double slots_used=((MQF*)dataStructure)->calculate_slotsUsedInCounting();
  //     slots_used-=((MQF*)dataStructure)->get_MQF()->metadata->noccupied_slots;
  //     results_slotsUsed.push_back(make_pair(capacity,slots_used));
  //   }
  //
  //
  // }
  // cout<<"Capacity\t1M Query per second"<<endl;
  // for(int i=0;i<result.size();i++)
  // {
  //   cout<<result[i].first<<"\t"<<result[i].second<<endl;
  // }
  // if(dataStrucureInput=="mqf"){
  //   cout<<"Capacity\t%slots used in counting"<<endl;
  //   for(int i=0;i<results_slotsUsed.size();i++)
  //   {
  //     cout<<results_slotsUsed[i].first<<"\t"<<results_slotsUsed[i].second<<endl;
  //   }
  // }
  // cerr<<"#generated elements= "<<g->get_generated_elements()<<endl;
  return 0;
}
