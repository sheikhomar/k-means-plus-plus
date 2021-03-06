#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <string>
#include <boost/array.hpp>
#include <blaze/Math.h>
#include <boost/range/algorithm_ext/erase.hpp>

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

void runParseEnronData() {
   Timer clock;
   clock.start();
   blaze::CompressedMatrix<uint> data = readEnronDataSet("data/docword.enron.txt");
   clock.stop();

   size_t durationInMs = clock.durationInMs();
   printf("Read %lu rows!\n", data.rows());
   printf("Runtime %lu ms (%lu secs)\n", durationInMs, durationInMs / 1000);
}

void testParseLineAsIntegerArray() {
   uint d[3] = {0, 0, 0};

   parseLineAsIntegerArray("4 1153 13\n", " ", d, 3);

   for (uint i = 0; i < 3; i++)
   {
      printf("%d\n", d[i]);
   }
}

blaze::DynamicMatrix<double> initialiseCentroids(const blaze::DynamicMatrix<double>& matrix, uint k) {
   size_t n = matrix.rows();
   size_t d = matrix.columns();

   // Initialise the sequence of pseudo-random numbers with a fixed random seed.
   static std::random_device seed;
   static std::mt19937 randomEngine(42);
   uniform_int_distribution<int> randomGen(0, n-1);

   blaze::DynamicMatrix<double> centrioids(k, d);

   for (size_t i = 0; i < k; i++) {
      int centroidIndex = randomGen(randomEngine);
      blaze::row(centrioids, i) = blaze::row(matrix, centroidIndex);
   }

   return centrioids;
}

blaze::DynamicMatrix<double> initialiseCentroidsKMeansPlusPlus(const blaze::DynamicMatrix<double>& matrix, uint k) {
   size_t n = matrix.rows();
   size_t d = matrix.columns();

   // Initialise the sequence of pseudo-random numbers with a fixed random seed.
   static std::random_device seed;
   static std::mt19937 randomEngine(42);
   uniform_int_distribution<int> randomGen(0, n-1);

   blaze::DynamicMatrix<double> centrioids(k, d);
   std::vector<uint> availableIndices(n);

   // Fill with 0, 1, 2, ..., N.
   std::iota(availableIndices.begin(), availableIndices.end(), 0); 
   
   uint chosenCentroids = 1;
   for (size_t c = 0; c < k; c++) {
      uint centroidIndex = -1;

      if (c == 0) {
         // Pick the first centroid uniformly at random.
         centroidIndex = randomGen(randomEngine);
      } else {

         blaze::DynamicVector<double> smallestDistances(n);

         // For each point, find the distance to the nearest centroid for all
         // the centroids that are select so far.
         for (uint p : availableIndices) {
            double smallestDistance = numeric_limits<double>::max();
            
            // Loop through previously selected clusters.
            for (size_t c2 = 0; c2 < c; c2++) {

               // Compute the L2 norm between point p and centroid c2.
               const double distance = blaze::norm(blaze::row(matrix, p) - blaze::row(centrioids, c2));

               // Decide if current distance is better.
               if (distance < smallestDistance) {
                  smallestDistance = distance;
               }
            }

            smallestDistances[p] = smallestDistance;
         }
         
         // Pick a point based on a weighted probability

         // Square distances
         smallestDistances *= smallestDistances;

         // Normalise.
         smallestDistances /= blaze::sum(smallestDistances);

         // Pick the index of a point randomly selected based on the weights.
         std::discrete_distribution<uint> weightedChoice(smallestDistances.begin(), smallestDistances.end());
         uint nextClusterCandidate = weightedChoice(randomEngine);

         // Assign centroid index.
         centroidIndex = availableIndices[nextClusterCandidate];
      }
      
      cout << "Centroid index for " << c << " => " << centroidIndex << "\n";
      
      // Copy point over centroids matrix.
      blaze::row(centrioids, c) = blaze::row(matrix, centroidIndex);

      // Remove it from the candidate list so the point cannot be
      // picked as another centroid.
      boost::remove_erase(availableIndices, centroidIndex);
   }

   return centrioids;
}

void kMeansLloyd(const blaze::DynamicMatrix<double>& matrix, blaze::DynamicMatrix<double>& centrioids, uint k, uint maxIterations) {
   size_t n = matrix.rows();
   size_t d = matrix.columns();

   // Initialise the sequence of pseudo-random numbers with a fixed random seed.
   static std::random_device seed;
   static std::mt19937 randomEngine(42);
   uniform_int_distribution<int> randomGen(0, n-1);

   blaze::DynamicVector<size_t> clusterAssignments(n);
   blaze::DynamicVector<size_t> clusterMemberCounts(k);

   for (size_t i = 0; i < maxIterations; i++) {
      // For each data point, assign the centroid that is closest to it.
      for (size_t p = 0; p < n; p++) {
         double bestDistance = numeric_limits<double>::max();
         size_t bestCluster = 0;
         
         // Loop through all the clusters.
         for (size_t c = 0; c < k; c++) {

            // Compute the L2 norm between point p and centroid c.
            const double distance = blaze::norm(blaze::row(matrix, p) - blaze::row(centrioids, c));

            // Decide if current distance is better.
            if (distance < bestDistance) {
               bestDistance = distance;
               bestCluster = c;
            }
         }

         // Assign cluster to point p.
         clusterAssignments[p] = bestCluster;
      }

      // Move centroids based on the cluster assignments.
      
      // First, save a copy of the centroids matrix.
      blaze::DynamicMatrix<double> oldCentrioids(centrioids);

      // Set all elements to zero.
      centrioids = 0; // Reset centroids.
      clusterMemberCounts = 0; // Reset cluster member counts.
      
      for (size_t p = 0; p < n; p++) {
         const size_t c = clusterAssignments[p];
         blaze::row(centrioids, c) += blaze::row(matrix, p);
         clusterMemberCounts[c] += 1;
      }

      for (size_t c = 0; c < k; c++) {
         const auto count = std::max<size_t>(1, clusterMemberCounts[c]);
         blaze::row(centrioids, c) /= count;
      }

      cout << "Centroids after iteration " << i << ": \n" << centrioids << "\n";

      auto diff = blaze::sum(blaze::abs(centrioids - oldCentrioids));

      if (diff < 0.001) {
         cout << "Stopping k-Means as centroids do not improve: " << diff << "\n";
         break;
      }
   }
}

void runKMeansLloyd() {
   blaze::DynamicMatrix<double> data{
      { -0.794152276623841F, 2.104951171962879F, },
      { -9.151551856068068F, -4.812864488195191F, },
      { -11.44182630908269F, -4.4578144103096555F, },
      { -9.767617768288718F, -3.19133736705118F, },
      { -4.536556476851341F, -8.401862882339504F, },
      { -6.263021151786394F, -8.1066608061999F, },
      { -6.384812343779634F, -8.473029703522716F, },
      { -9.204905637733754F, -4.5768792770429965F, },
      { -2.760179083161441F, 5.551213578682775F, },
      { -1.1710417594110225F, 4.330918155822106F, },
      { -10.036408012919152F, -5.5691209020665F, },
      { -9.875891232661665F, -2.823864639451285F, },
      { -7.175329210075055F, -8.770590168336406F, },
      { -2.406718199699357F, 6.098944469870908F, },
      { -4.874182454688006F, -10.049589027515138F, },
      { -6.078546995836497F, -7.939694203288603F, },
      { -6.832387624479001F, -7.47067669775956F, },
      { -2.346732606068119F, 3.561284227344442F, },
      { -10.341566179224177F, -3.9097516905289575F, },
      { -11.092624349394143F, -3.7839661143045364F, },
      { -6.502121087038712F, -7.912491012386313F, },
      { -10.263931009656723F, -3.920734000669846F, },
      { -6.816083022269968F, -8.449869256994909F, },
      { -1.340520809891421F, 4.157119493365752F, },
      { -10.372997453743215F, -4.592078954817427F, },
      { -7.374998957175799F, -10.588065868731183F, },
      { -6.623517740089157F, -8.253383337907545F, },
      { -1.359389585992692F, 4.054240022349643F, },
      { -0.19745196890354544F, 2.3463491593455075F, },
      { -6.5443058465843675F, -9.297569494247188F, },
      { -1.9274479855745354F, 4.9368453355813475F, },
      { -2.8020781039706595F, 4.057147146430284F, },
      { -7.581976641577967F, -9.150254932274308F, },
      { -1.8513954583101344F, 3.5188609047583252F, },
      { -8.370061750504195F, -3.615336850788729F, },
      { -7.251451955565088F, -8.25497397715319F, },
      { -8.798794623751593F, -3.7681921298792607F, },
      { -11.370829823899857F, -3.6381891553209127F, },
      { -10.178632805731251F, -4.557269175156462F, },
      { -7.2013269275537715F, -8.272282292398854F, },
      { -6.7842171065351F, -8.226340808371322F, },
      { -9.647166524988995F, -5.265631958600636F, },
      { -1.9819771099620271F, 4.022435514174746F, },
      { -11.227770639320063F, -3.402811051386989F, },
      { -9.799412783526332F, -3.834339901555746F, },
      { -6.5354168593050295F, -8.015526894626658F, },
      { -0.757969185355724F, 4.908984207745029F, },
      { 0.5260155005846419F, 3.009993533355024F, },
      { -2.7768702545837973F, 4.640905566660254F, },
      { -1.7824501314671677F, 3.4707204345840927F, },
      { -10.220040646263461F, -4.154106616293202F, },
      { -6.4058323875575285F, -9.780666445240302F, },
      { -6.987061055745032F, -7.5348478426255205F, },
      { -7.465760375446665F, -7.329222486173637F, },
      { -1.5394009534668904F, 5.023692978550581F, },
      { -6.569670859679778F, -8.327931264366546F, },
      { -10.617713347601232F, -3.255316513290986F, },
      { -8.723956573494325F, -1.98624679810847F, },
      { -1.6173461592329268F, 4.9893050825589835F, },
      { -1.1466300855305107F, 4.108397033740446F, },
      { -9.811151112664817F, -3.543296900154948F, },
      { -7.711798871912647F, -7.251741212975334F, },
      { -6.561697370222412F, -6.860002222091783F, },
      { -10.02232945952888F, -4.728510166532364F, },
      { -11.855694368099854F, -2.7171845169103843F, },
      { -5.733425071070147F, -8.440535968100065F, },
      { -2.4139578469451726F, 5.659358024076449F, },
      { -8.337440938903733F, -7.839680384160613F, },
      { -1.8319881134989553F, 3.5286314509217895F, },
      { -9.574218149588988F, -3.8760084790146454F, },
      { -9.5942208618623F, -3.3597700241261377F, },
      { -9.257156052556827F, -4.907049149171139F, },
      { -6.46256290336211F, -7.732945900976985F, },
      { -0.8205764920740146F, 5.337591950146718F, },
      { 0.00024227116135100424F, 5.148534029420497F, },
      { -9.682077556411496F, -5.975549763187208F, },
      { -6.195996026651871F, -7.402816464759037F, },
      { -7.02121319047935F, -8.379542347137651F, },
      { -2.187731658211975F, 3.333521246686991F, },
      { -10.4448410684391F, -2.7288408425577058F, },
      { -0.5279305184970926F, 5.92630668526536F, },
      { -11.196980535988288F, -3.090003229819183F, },
      { -9.837675434205272F, -3.0771796262469797F, },
      { -5.160223475316758F, -7.04217140606354F, },
      { -2.351220657673829F, 4.0097363419871845F, },
      { -0.5257904636130821F, 3.3065986015291307F, },
      { -1.4686444212810534F, 6.506745005322004F, },
      { -0.7587039566841077F, 3.7227620096688283F, },
      { -10.303916516281474F, -3.1253739047559583F, },
      { -2.3308060367853387F, 4.39382526992426F, },
      { -5.904543613663969F, -7.783735388248322F, },
      { -1.6087521511724905F, 3.769494222273808F, },
      { -1.8684541393232976F, 4.993113060025359F, },
      { -10.668374789942131F, -3.5757847610422853F, },
      { -8.876294795417436F, -3.544448009426377F, },
      { -6.026057581798325F, -5.966248457649787F, },
      { -7.047472775357734F, -9.275246833370932F, },
      { -1.3739725806942609F, 5.291631033113889F, },
      { -6.2539305108541825F, -7.108786009916786F, },
      { 0.08525185826796045F, 3.6452829679480585F, },
   };
   
   auto centrioids = initialiseCentroidsKMeansPlusPlus(data, 3);
   //auto centrioids = initialiseCentroids(data, 3);
   cout << "Initial centroids: \n" << centrioids << "\n";

   kMeansLloyd(data, centrioids, 3, 20);
}

int main()
{
   // runParseEnronData();
   runKMeansLloyd();
   cout << "Done!\n";
}
