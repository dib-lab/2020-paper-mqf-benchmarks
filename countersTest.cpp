#include "gqf.h"
#include "cqf/gqf.h"
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

    noSlots = (uint64_t) ((double) noSlots * 0.90);
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
                    *res_noSlots = noSlots;
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


  uint64_t p=8;
  uint64_t qbits=4;
  uint64_t slot_size=4;
  uint64_t fixedCounterSize=2;
  QF mqf;
  qf_init(&mqf, (1ULL<<qbits), qbits+slot_size, 0,fixedCounterSize, 0,true, "", 2038074761);

  uint64_t c=1;
  int xaxis=16;
  int yaxis=15;

  cout<<"Hash reminder\tCount\tLog(MQF bits/CQF bits)\n";
  for(int j=0;j<yaxis;j++) {
      for (int i = 0; i < xaxis; i++) {
          uint64_t mqfSlots = mqf.metadata->noccupied_slots;
          qf_insert(&mqf, i, c);
          int mqfBits=(mqf.metadata->noccupied_slots - mqfSlots)*6;
          int mqfSolts=(mqf.metadata->noccupied_slots - mqfSlots);
          //cout << "\t" << (mqf.metadata->noccupied_slots - mqfSlots)*6;
          qf_remove(&mqf, i, c);


          cqf::QF ccqf;
          cqf::qf_init(&ccqf, (1ULL << qbits), qbits + slot_size, 0, true, "", 2038074761);


          cqf::qf_insert(&ccqf, i, 0, c, false, false);
          //cout << "\t" << ccqf.metadata->noccupied_slots*4 << endl;
          int cqfBits=ccqf.metadata->noccupied_slots*4;

          cout << i << "\t" << c<<"\t"<<log((double)mqfBits/(double)cqfBits)<<"\t"<<mqfSolts<<"\n";
          cqf::qf_destroy(&ccqf, true);

      }
      c*=2;
      //cout<<endl;
  }


    return 0;
}
