#include "gqf.h"
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>
#include<iostream>
#include <fstream>
#include "hashutil.h"
#include "kmer.h"
using namespace std;

int main(int argc, char const *argv[]) {

  ifstream dataset(argv[1]);
  uint64_t qbits=atoi(argv[2]);
  uint64_t nrepeats=atoi(argv[3]);
  QF qf;
  srand (1);

  uint64_t num_hash_bits=qbits+8;

  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,4, true, "", 2038074761);
  string kmer;
  uint64_t countedKmers=0;
  while(dataset>>kmer)
  {
    uint64_t item=kmercounting::str_to_int(kmer);
    uint64_t hash=kmercounting::HashUtil::MurmurHash64A(((void*)&item), sizeof(item),qf.metadata->seed);
    hash=hash%qf.metadata->range;


    if(qf_count_key(&qf, hash)==0){
      qf_insert(&qf,hash,nrepeats,false,false);
      countedKmers++;
      cout<<nrepeats<<"\t"<<countedKmers<<endl;
    }
  }
//  qf_destroy(&qf,true);
  return 0;
}
