//
//    AEM Invert : Software for inversion of AEM data using the
//    trans-dimensional tree method and forward modelling code
//    written by Ross Brodie from Geoscience Australia. See
//
//      R Hawkins, R Brodie and M Sambridge, "Bayesian trans-dimensional inversion of
//    Airborne Electromagnetic 2D Conductivity profiles", Exploration Geophysics, 2017
//    https://doi.org/10.1071/EG16139
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <getopt.h>

extern "C" {
#include "chain_history.h"
  
};

#include "global.hpp"

#include "aemutil.hpp"

#include <map>

struct coefficient_counter {

  coefficient_counter() :
    pv(0),
    av(0),
    pb(0),
    ab(0),
    pd(0),
    ad(0),
    meann(0),
    mean(0.0),
    var(0.0)
  {
  }
  
  int pv;
  int av;

  int pb;
  int ab;

  int pd;
  int ad;

  int meann;
  double mean;
  double var;
};

struct user_data {
  int thincounter;
  int thin;
  int skip;
  
  int counter;
  
  std::map<int, coefficient_counter*> coefficients;
};

static int process(int i,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v);

static char short_options[] = "i:o:t:s:S:h";
static struct option long_options[] = {

  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},

  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"maxsteps", required_argument, 0, 'S'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  chain_history_t *ch;
  
  std::vector<std::string> input_file;
  
  char *output_file;

  int thin;
  int skip;
  int maxsteps;

  FILE *fp_in;
  FILE *fp_out;

  struct user_data data;
  multiset_int_double_t *S_v;

  int i;
  int j;

  int processesperchain;
  
  int mpi_size;
  int mpi_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    

  /*
   * Default values
   */
  fp_in = NULL;
  fp_out = NULL;
  
  output_file = NULL;

  thin = 0;
  skip = 0;

  maxsteps = 1000000;

  processesperchain = 1;
  
  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {
    case 'i':
      input_file.push_back(optarg);
      break;

    case 'o':
      output_file = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      break;

    case 's':
      skip = atoi(optarg);
      break;
      
    case 'S':
      maxsteps = atoi(optarg);
      if (maxsteps < 1000) {
	fprintf(stderr, "error: maxsteps should be 1000 or greater\n");
	return -1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
      
    }
  }

  if (input_file.size() == 0) {
    fprintf(stderr, "error: required parameter input file missing\n");
    return -1;
  }

  if (output_file == NULL) {
    fprintf(stderr, "error: required parameter output file missing\n");
    return -1;
  }

  ch = chain_history_create(maxsteps);
  if (ch == NULL) {
    fprintf(stderr, "error: failed to create chain history\n");
    return -1;
  }
  
  data.thincounter = 0;
  data.thin = thin;
  data.skip = skip;

  data.counter = 0;

  S_v = multiset_int_double_create();
  if (S_v == NULL) {
    fprintf(stderr, "error: failed to create multiset\n");
    return -1;
  }

  for (auto &infile: input_file) {

    std::string chfile = mkfilenamerank(nullptr, infile.c_str(), mpi_rank * processesperchain);
    fp_in = fopen(chfile.c_str(), "r");
    if (fp_in == NULL) {
      fprintf(stderr, "error: failed to open input file\n");
      return -1;
    }
    printf("Loaded: %s\n", chfile.c_str());
    
    /*
     * Process the chain history
     */
    while (!feof(fp_in)) {
      
      if (chain_history_read(ch,
			     (ch_read_t)fread,
			     fp_in) < 0) {
	if (feof(fp_in)) {
	  break;
	}
	
	fprintf(stderr, "error: failed to read chain history\n");
	return -1;
      }
      
      if (chain_history_replay(ch,
			       S_v,
			       (chain_history_replay_function_t)process,
			       &data) < 0) {
	fprintf(stderr, "error: failed to replay\n");
	return -1;
      }
    }
    printf("%d records\n", data.counter);
    fclose(fp_in);
  }
    
  MPI_Barrier(MPI_COMM_WORLD);
  
  chain_history_destroy(ch);
  multiset_int_double_destroy(S_v);

  delete [] data.mean;
  delete [] data.variance;
  delete [] data.model;
  delete [] data.workspace;

  MPI_Finalize();
  
  return 0;
}

static int process(int stepi,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v)
{
  struct user_data *d = (struct user_data *)user;
  double delta;
  int i;
  int hi;
  
  if ((d->thincounter >= d->skip) && (d->thin <= 1 || (d->thincounter % d->thin) == 0)) {

    
  }
  d->thincounter ++;
  
  return 0;
}


static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <file>                Input ch file\n"
	  " -o|--output <file>               Output mean model file\n"
	  "\n"
	  " -t|--thin <int>                  Only processing every ith sample\n"
	  " -s|--skip <int>                  Skip n samples from beginning\n"
	  "\n"
	  " -S|--maxsteps <int>              Chain history max steps\n"
	  "\n"
	  " -h|--help            Show usage\n"
	  "\n",
	  pname);
}

 
