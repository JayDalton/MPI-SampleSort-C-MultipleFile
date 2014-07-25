#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#define LSIZE 24000000  // Lines of File
#define TSIZE 32  // Structur size
//#define FOUT "/home/jay/Schreibtisch/SampleSort/DataSmall/twitter.out"
//#define FLOG "/home/jay/Schreibtisch/SampleSort/DataSmall/twitter.log"
#define FOUT "/mpidata/ergebnisse/rhtg_twitter.out"
#define FLOG "/mpidata/ergebnisse/rhtg_twitter.log"

FILE* LogFile;
char* MONTHS[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
char* FILES[] = {
  "/mpidata/parsys14/gross/twitter.data.0", 
  "/mpidata/parsys14/gross/twitter.data.1", 
  "/mpidata/parsys14/gross/twitter.data.2", 
  "/mpidata/parsys14/gross/twitter.data.3", 
  "/mpidata/parsys14/gross/twitter.data.4", 
  "/mpidata/parsys14/gross/twitter.data.5", 
  "/mpidata/parsys14/gross/twitter.data.6", 
  "/mpidata/parsys14/gross/twitter.data.7", 
  "/mpidata/parsys14/gross/twitter.data.8", 
  "/mpidata/parsys14/gross/twitter.data.9", 
  "/mpidata/parsys14/gross/twitter.data.10", 
  "/mpidata/parsys14/gross/twitter.data.11", 
  "/mpidata/parsys14/gross/twitter.data.12", 
  "/mpidata/parsys14/gross/twitter.data.13", 
  "/mpidata/parsys14/gross/twitter.data.14", 
  "/mpidata/parsys14/gross/twitter.data.15"
};

int readNumber(char** lptr) {
  char* ptr = *lptr;
  char* line = *lptr;
  while (*ptr != ' ') {
    ptr++;
  }
  *ptr = 0;
  *lptr = ptr + 1;
  return atoi(line);
}

int readMonth(char** lptr) {
  char* ptr = *lptr;
  char* line = *lptr;
  while (*ptr != ' ') {
    ptr++;
  }
  *ptr = 0;
  *lptr = ptr + 1;
  int i, m;
  for (i = 0, m = 1; i < 12; i++, m++) {
    if (strncmp(line, MONTHS[i], 3) == 0) {
      return m;
    }
  }
  fprintf(stderr, "invalid month: %s\n", line);
  exit(3);
}

int countHits(const char* line, const char* key) 
{
  int n = strlen(key);
  int k = strlen(line) - n;
  int i;
  int hits = 0;
  for (i = 0; i < k; i++, line++) {
    if (*line == *key) {
      if (strncmp(line, key, n) == 0) {
        hits++;
      }
    }
  }
  return hits;
}

union {
  unsigned long longint;
  unsigned char byte[8];
} LongConv;

void printTweet(const char* t) 
{
  fprintf(LogFile, "Tweet[");
  int i;
  for (i = 0; i < TSIZE; i++) {
    int k = t[i];
    fprintf(LogFile, "%2x ", k < 0 ? k + 256 : k);
  }
  fprintf(LogFile, "]\n");
}

void printTweetList(char* liste, unsigned long anzahl)
{
  unsigned long i;
//  for( i = 0; i < anzahl; i++) {
//    printTweet(liste + (i * TSIZE));
//  }
}

void writeTweet(char* tweet, const int fn, const int ln, const int hits,
    const int month, const int day, char* line) {
  short* ptr1 = (short*) tweet;
  *ptr1 = (short) fn;
  int* ptr2 = (int*) (tweet + 2);
  *ptr2 = ln;
  *(tweet + 6) = (char) hits;
  *(tweet + 7) = (char) month;
  *(tweet + 8) = (char) day;
  int i;
  int n = TSIZE - 9;
  for (i = strlen(line); i < n; i++)
    line[i] = ' '; // padding
  memcpy(tweet + 9, line, n);
       //printTweet(tweet); printf("\n");
}

void saveTweet(char* tweet, const int fn, const int ln, const int hits, const int month, const int day, char* line) 
{
  short* ptr1 = (short*) tweet;
  *ptr1 = (short) fn;
  int* ptr2 = (int*) (tweet + 2);
  *ptr2 = ln;
  *(tweet + 6) = (char) hits;
  *(tweet + 7) = (char) month;
  *(tweet + 8) = (char) day;
  int i;
  int n = TSIZE - 9;
  for (i = strlen(line); i < n; i++) {
    line[i] = ' '; // padding
  }
  memcpy(tweet + 9, line, n);
}

void writeTweets(char* output, const char* fnam, const unsigned long lines) 
{
	FILE* file = fopen(fnam, "w");
  if (file == NULL) {
    fprintf(stderr, "open failed: %s\n", fnam);
  }
	unsigned long i;
	char* tweet;
	for (i = 0, tweet = output; i < lines; i++, tweet += TSIZE) {
		short* fnp = (short*) tweet;
		int* lnp = (int*) (tweet + 2);
		fprintf(file, "%d %d\n", *fnp, *lnp);
	}
	fclose(file);
}

void loadTweets(char* input, const char* fnam, const unsigned long lines, const char* key) 
{
  FILE* file = fopen(fnam, "r");
  if (file == NULL) {
    fprintf(stderr, "open failed: %s\n", fnam);
  }

  unsigned long i;
  char buffer[1024];
  char* tweet;

  for (i = 0, tweet = input; i < lines; i++, tweet += TSIZE) {
    char* line = fgets(buffer, 1024, file);
    if (line == NULL) {
      fprintf(stderr, "error reading line %lu\n", i);
      exit(2);
    }
    int fn = readNumber(&line);
    int ln = readNumber(&line);
    int month = readMonth(&line);
    int day = readNumber(&line);
    int hits = countHits(line, key);
    saveTweet(tweet, fn, ln, hits, month, day, line);
  }
  fclose(file);
}

int compTweets(const void* t1, const void* t2)
{
  return memcmp(t2 + 6, t1 + 6, TSIZE - 6);
}

void copyTweet(char *from, char *to) 
{
  memcpy(to, from, TSIZE * sizeof(unsigned char));
}

void copyTweetList(char *from, char *to, unsigned long anzahl)
{
  memcpy(to, from, anzahl * TSIZE * sizeof(unsigned char));
}

/**** Function Declaration Section ****/

main (int argc, char *argv[])
{
  /* Variable Declarations */

  int             Numprocs, MyRank, Root = 0;
  unsigned long   i, j, k, count, NoElementsToSort;
  unsigned char   *InputData, *Tweets;
  unsigned char   FileName[1024];
  unsigned char   *Splitter, *AllSplitter;
  unsigned char   *Buckets, *BucketBuffer, *LocalBucket;
  unsigned char   *OutputBuffer, *Output;
  FILE            *InputFile, *fp;
  MPI_Status      status; 
  
  /**** Initialising ****/
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

  double starttime = MPI_Wtime();

  unsigned long FilesSizeAll = sizeof(FILES) / sizeof(FILES[0]);  // Anzahl aller Dateien
  unsigned long FilesSizeNode = 1;//FilesSizeAll / Numprocs;          // Anzahl Dateien pro Prozess
  unsigned long TweetSizeFile = LSIZE;                            // Anzahl Zeilen pro Datei
  unsigned long TweetSizeNode = TweetSizeFile;// * FilesSizeNode;    // Anzahl Tweets pro Node
  unsigned long TweetSizeAll = TweetSizeNode * Numprocs;          // Anzahl Tweets insgesamt
    
  if (argc != 2) {
    if (MyRank == 0) { 
      printf(" Usage : run size\n");
    }
    MPI_Finalize();
    exit(0);
  }

  if (( TweetSizeAll % Numprocs) != 0){
    if (MyRank == Root) {
      fprintf(LogFile, "Sending Data: Number of Elements are not divisible by Numprocs \r\n");
    }
    MPI_Finalize();
    exit(0);
  }


  // Init Logging
  sprintf(FileName, "%s.%s.%d", FLOG, argv[1], MyRank);
  LogFile = fopen(FileName, "w");
  
  fprintf(LogFile, "********* Statistics *********\n");
  fprintf(LogFile, "Numprocs: %d\n", Numprocs);
  fprintf(LogFile, "MyRank: %d\n", MyRank);
  fprintf(LogFile, "FilesSizeAll: %lu\n", FilesSizeAll);
  fprintf(LogFile, "FilesSizeNode: %lu\n", FilesSizeNode);
  fprintf(LogFile, "TweetSizeFile: %lu\n", TweetSizeFile);
  fprintf(LogFile, "TweetSizeNode: %lu\n", TweetSizeNode);
  fprintf(LogFile, "TweetSizeNode: %lu Bytes\n", TweetSizeNode * TSIZE);
  fprintf(LogFile, "TweetSizeAll: %lu\n", TweetSizeAll);
  fprintf(LogFile, "TweetSizeAll: %lu Bytes\n", TweetSizeAll * TSIZE);
  fprintf(LogFile, "Time Initial: %f\n", MPI_Wtime() - starttime);
  fprintf(LogFile, "********* Statistics *********\n\n");
  starttime = MPI_Wtime();

  /**** Reading Input Parallel ****/
  fprintf(LogFile, "Each Node reading local Input: %lu Bytes\n", TweetSizeNode * TSIZE);
  InputData = (unsigned char*) calloc (TweetSizeNode, TSIZE * sizeof(unsigned char));

  if (InputData == NULL) {
    fprintf(LogFile, "Error : Can not allocate memory \n");
  }

  loadTweets(InputData, FILES[MyRank], TweetSizeFile, argv[1]);
  fprintf(LogFile, "File: %lu ID: %i Name: %s\n", i, MyRank, FILES[MyRank]);
  
  fprintf(LogFile, "Local Tweet List: %lu\n", TweetSizeNode);
  printTweetList(InputData, TweetSizeNode);

  fprintf(LogFile, "Time Load Input: %f\n", MPI_Wtime() - starttime);
  starttime = MPI_Wtime();
  
  /**** Sorting Locally ****/
  qsort ((unsigned char *)InputData, TweetSizeNode, TSIZE, compTweets);
  fprintf(LogFile, "InputData sorted: %lu\n", TweetSizeNode);
  printTweetList(InputData, TweetSizeNode);

  /**** Choosing Local Splitters ****/
  unsigned long NumSplitter = Numprocs - 1;
  unsigned long TweetBucketSize = TweetSizeNode / Numprocs;

  fprintf(LogFile, "* Choosing Local Splitters: %lu, TweetBucketSize: %lu  * \n", NumSplitter, TweetBucketSize);
  Splitter = (unsigned char *) calloc (NumSplitter, TSIZE * sizeof(unsigned char));

  for (i = 0; i < NumSplitter; i++) 
  {
    copyTweet(InputData + (i + 1) * (TweetBucketSize * TSIZE), Splitter + (i * TSIZE));
  } 

  fprintf(LogFile, "Splitter local: %lu \n",  NumSplitter);
  printTweetList(Splitter, NumSplitter);

  /**** Gathering Local Splitters at Root ****/
  fprintf(LogFile, "* Gathering Local Splitters at Root * \n");
  AllSplitter = (unsigned char*) calloc (Numprocs * NumSplitter, TSIZE * sizeof (unsigned char));

  MPI_Gather (
    Splitter, NumSplitter * TSIZE, MPI_UNSIGNED_CHAR, 
    AllSplitter, NumSplitter * TSIZE, MPI_UNSIGNED_CHAR, 
    Root, MPI_COMM_WORLD
  );

  fprintf(LogFile, "AllSplitter: local %lu\n", Numprocs * NumSplitter);
  printTweetList(AllSplitter, Numprocs * NumSplitter);

  /**** Choosing Global Splitters ****/
  if (MyRank == Root) {

    qsort ((unsigned char *) AllSplitter, Numprocs * NumSplitter, TSIZE, compTweets);

    fprintf(LogFile, "All Splitter Root: sorted %lu\n", Numprocs * NumSplitter);
    printTweetList(AllSplitter, Numprocs * NumSplitter);

    for (i = 0; i < NumSplitter; i++) 
    {
      copyTweet(AllSplitter + (i + 1) * (NumSplitter * TSIZE), Splitter + (i * TSIZE));
    } 

    fprintf(LogFile, "Splitter Root selected: %lu\n", NumSplitter);
    printTweetList(Splitter, NumSplitter);
  }
  
  /**** Broadcasting Global Splitters ****/
  fprintf(LogFile, "Broadcasting Global Splitters\n");
  MPI_Bcast (Splitter, NumSplitter * TSIZE, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  fprintf(LogFile, "Splitter local: %lu\n", NumSplitter);
  printTweetList(Splitter, NumSplitter);

  /**** Creating Numprocs Buckets locally ****/
  fprintf(LogFile, "Creating Numprocs Buckets locally: %lu + %i\n", TweetSizeNode, Numprocs);
  Buckets = (unsigned char *) calloc (TweetSizeAll + Numprocs, TSIZE * sizeof(unsigned char));

  i = 0;  // Index InputData
  j = 0;  // Index Splitter/Bucket
  k = 1;  // Anzahl Tweets pro Bucket

  for (i = 0; i < TweetSizeNode; i++) {
    
    if (j < NumSplitter) {
      
      if (0 > compTweets(InputData + (i * TSIZE), Splitter + (j * TSIZE))) {

        copyTweet(InputData + (i * TSIZE), Buckets + (((TweetSizeNode + 1) * j) + k) * TSIZE);
        k++;
      
      } else {
        LongConv.longint = k - 1;
        fprintf(LogFile, "Anzahl: %lu ", LongConv.longint);
        memcpy(Buckets + ((TweetSizeNode + 1) * j * TSIZE), LongConv.byte, sizeof(LongConv));
		    k=1;
			  j++;
		    i--;
       }
    } else {
      copyTweet(InputData + (i * TSIZE), Buckets + (((TweetSizeNode + 1) * j) + k) * TSIZE);
      k++;
	  }
  }
  LongConv.longint = k - 1;
  fprintf(LogFile, "Anzahl: %lu\n", LongConv.longint);
  memcpy(Buckets + ((TweetSizeNode + 1) * j * TSIZE), LongConv.byte, sizeof(LongConv));

  fprintf(LogFile, "Local Buckets List: %lu + %i\n", TweetSizeNode, Numprocs);
  printTweetList(Buckets, TweetSizeAll + Numprocs);
  
  free(InputData);

  /**** Sending buckets to respective processors ****/
  fprintf(LogFile, "Sending buckets to respective processors\n");
  BucketBuffer = (unsigned char *) calloc (TweetSizeAll + Numprocs, TSIZE * sizeof (unsigned char));

  MPI_Alltoall (  // MPI_Allgather is same?
    Buckets, (TweetSizeNode + 1) * TSIZE, MPI_UNSIGNED_CHAR, 
    BucketBuffer, (TweetSizeNode + 1) * TSIZE, MPI_UNSIGNED_CHAR, 
    MPI_COMM_WORLD
  );
    
  free(Buckets);
  
  fprintf(LogFile, "BucketBufferAll List: %lu\n", TweetSizeAll + Numprocs);
  printTweetList(BucketBuffer, TweetSizeAll + Numprocs);

  /**** Rearranging BucketBuffer ****/
  fprintf(LogFile, "Rearranging BucketBuffer\n");
  LocalBucket = (unsigned char *) calloc (2 * TweetSizeAll / Numprocs,  TSIZE * sizeof(unsigned char));
  fprintf(LogFile, "LocalBucket Alloc: %lu\n", 2 * TweetSizeAll / Numprocs);

  count = 1;

  for (j = 0; j < Numprocs; j++) 
  {
    k = 1;
    memcpy(LongConv.byte, BucketBuffer + (TweetSizeNode + 1) * j * TSIZE, sizeof(LongConv));
    fprintf(LogFile, "Anzahl: %lu ", LongConv.longint);
    for (i = 0; i < LongConv.longint; i++) 
    {
      copyTweet(BucketBuffer + (((TweetSizeNode + 1) * j + k)) * TSIZE, LocalBucket + (count * TSIZE));
      count++;
      k++;
    }
  }
  LongConv.longint = count - 1;
  fprintf(LogFile, "Anzahl: %lu\n", LongConv.longint);
  memcpy(LocalBucket, LongConv.byte, sizeof(LongConv));

  fprintf(LogFile, "LocalBucket List: %lu\n", 2 * TweetSizeAll);
  printTweetList(LocalBucket, 2 * TweetSizeAll);

  free(BucketBuffer);

  /**** Sorting Local Buckets using Bubble Sort ****/
  /*qsort ((char *) InputData, NoofElements_Bloc, sizeof(int), intcompare); */
  fprintf(LogFile, "Sorting Local Buckets using Bubble Sort\n");

  NoElementsToSort = LocalBucket[0];
  memcpy(LongConv.byte, LocalBucket, sizeof(LongConv));
  NoElementsToSort = LongConv.longint;
  fprintf(LogFile, "NoElementsToSort: %lu\n", NoElementsToSort);

  qsort ((unsigned char *) LocalBucket + TSIZE, NoElementsToSort, TSIZE, compTweets); 
  fprintf(LogFile, "LocalBucket List sorted: %lu\n", 2 * TweetBucketSize);

  fprintf(LogFile, "Time Sort Data: %f\n", MPI_Wtime() - starttime);
  starttime = MPI_Wtime();

    /**** Printng the output ****/
  fprintf(LogFile, "**** Printng the output ****\n");
  printTweetList(LocalBucket, 2 * TweetBucketSize);

	sprintf(FileName, "%s.%s.%d", FOUT, argv[1], MyRank);
  writeTweets(LocalBucket + TSIZE, FileName, LongConv.longint);

  fprintf(LogFile, "Time File Writing: %f\n", MPI_Wtime() - starttime);
  starttime = MPI_Wtime();

//  free(InputData);
  free(Splitter);
  free(AllSplitter);
//  free(Buckets);
//  free(BucketBuffer);
  free(LocalBucket);

  fclose(LogFile);

  /**** Finalize ****/
  MPI_Finalize();
}
