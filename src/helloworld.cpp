#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <boost/array.hpp>
#include <blaze/Math.h>
using blaze::StaticVector;
using blaze::DynamicVector;

using namespace std;

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

void readEnronDataSet(const std::string& path) {
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

      uint counter = 0;
      while (getline(fileHandle, line)) {
         printf("Line: %s\n", line.c_str());

         uint numbers[3] = { 0, 0, 0 };
         parseLineAsIntegerArray(line, " ", numbers, 3);
         
         printf("- DocId:  %d\n", numbers[0]);
         printf("- WordId: %d\n", numbers[1]);
         printf("- Count:  %d\n", numbers[2]);

         counter++;

         if (counter > 10)
            break;
      }
      fileHandle.close();
   } else {
      cout << "File never opened!\n";
   }
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

   readEnronDataSet("data/docword.enron.txt");

   uint d[3] = {0, 0, 0};

   parseLineAsIntegerArray("4 1153 13\n", " ", d, 3);

   for (uint i = 0; i < 3; i++)
   {
      printf("%d\n", d[i]);
   }
   

   cout << "Done!\n";
}
