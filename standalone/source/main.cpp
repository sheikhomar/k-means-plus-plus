#include <cxxopts.hpp>

using namespace std;

auto main(int argc, char** argv) -> int {
    cxxopts::Options options(*argv, "A program to welcome the world!");
    
    cout << "Hello World! \n" ;
}