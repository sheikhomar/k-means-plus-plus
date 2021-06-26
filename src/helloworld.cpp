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

void readEnronDataSet(const std::string& path) {
   cout << "Reading Enron dataset from " << path << "...\n";
   
   ifstream fileHandle;
   fileHandle.open(path.c_str(), ios::in);
   if (fileHandle.is_open()) {
      std::string line;
      uint counter = 0;
      while (getline(fileHandle, line)) {
         printf("%s\n", line.c_str());
         counter++;
         if (counter > 10) {
            break;
         }
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
   cout << "Done!\n";
}
