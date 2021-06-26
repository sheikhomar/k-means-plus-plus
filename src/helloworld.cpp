#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <string>
#include <boost/array.hpp>
#include <blaze/Math.h>

using blaze::StaticVector;
using blaze::DynamicVector;

using namespace std;

template <class TimeT = std::chrono::milliseconds,
          class ClockT = std::chrono::steady_clock>
class Timer
{
    using timep_t = typename ClockT::time_point;
    timep_t _start = ClockT::now(), _end = {};

public:
    void start() { 
        _end = timep_t{}; 
        _start = ClockT::now(); 
    }
    
    void stop() { _end = ClockT::now(); }
    
    size_t durationInMs() const { 
        return (std::chrono::duration_cast<TimeT>(_end - _start).count());
    }
};



void parseLineAsIntegerArray(const std::string& line, const std::string& delimiter, uint* storage, uint storageSize) {
   int start = 0;
   int end = line.find(delimiter);
   int i = 0;
   while (end != -1) {
      int parsedInt = std::stoi(line.substr(start, end - start));
      if (i >= storageSize) {
         throw "Found more integers than allocated memory!";
      }
      storage[i] = parsedInt;
      i++;
      start = end + delimiter.size();
      end = line.find(delimiter, start);
   }

   int parsedInt = std::stoi(line.substr(start, end - start));
   if (i >= storageSize) {
      throw "Found more integers than allocated memory!";
   }
   storage[i] = parsedInt;
}

blaze::CompressedMatrix<uint> readEnronDataSet(const std::string& path) {
   cout << "Reading Enron dataset from " << path << "...\n";
   
   ifstream fileHandle;
   fileHandle.open(path.c_str(), ios::in);
   if (fileHandle.is_open()) {
      string line;

      // First line contains the number of documents.
      getline(fileHandle, line);
      uint dbSize = std::stoi(line);

      // Second line contains the number of words in the vocabulary
      getline(fileHandle, line);
      uint vocabularySize = std::stoi(line);

      // Third line contains the total number of words in the collection.
      // Ignore this for now.
      getline(fileHandle, line);

      printf("DB Size: %d\n", dbSize);
      printf("Vocab Size: %d\n", vocabularySize);

      blaze::CompressedMatrix<uint> X(dbSize, vocabularySize);

      uint counter = 0;
      while (getline(fileHandle, line)) {
         // printf("Line: %s\n", line.c_str());

         uint numbers[3] = { 0, 0, 0 };
         parseLineAsIntegerArray(line, " ", numbers, 3);

         uint docId = numbers[0];
         uint wordId = numbers[1];
         uint count = numbers[2];
         
         // printf("- DocId:  %d\n", docId);
         // printf("- WordId: %d\n", wordId);
         // printf("- Count:  %d\n", count);

         X(docId, wordId) = count;
         
         counter++;

         if (counter > 100000)
            break;
      }
      fileHandle.close();

      return X;
   } else {
      cout << "File never opened!\n";
   }
   blaze::CompressedMatrix<uint> zero(0, 0);
   return zero;
}

int main()
{
   vector<string> msg {"Hello", "World"};
   
   for (const string& word : msg)
   {
      cout << word << " ";
   }
   cout << endl;


   // Instantiation of a static 3D column vector. The vector is directly initialized as
   //    ( 4 -2  5 )
   StaticVector<int,3UL> a{ 4, -2, 5 };

   // Instantiation of a dynamic 3D column vector. Via the subscript operator the values are set to
   //    ( 2  5 -3 )
   DynamicVector<int> b( 3UL );
   b[0] = 2;
   b[1] = 5;
   b[2] = -3;

   // Adding the vectors a and b
   DynamicVector<int> c = a + b;

   // Printing the result of the vector addition
   std::cout << "c =\n" << c << "\n";

   Timer clock;
   clock.start();
   blaze::CompressedMatrix<uint> data = readEnronDataSet("data/docword.enron.txt");
   clock.stop();

   size_t durationInMs = clock.durationInMs();
   printf("Read %lu rows!\n", data.rows());
   printf("Runtime %lu ms (%lu secs)\n", durationInMs, durationInMs / 1000);

   uint d[3] = {0, 0, 0};

   parseLineAsIntegerArray("4 1153 13\n", " ", d, 3);

   for (uint i = 0; i < 3; i++)
   {
      printf("%d\n", d[i]);
   }
   

   cout << "Done!\n";
}
