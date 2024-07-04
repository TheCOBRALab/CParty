# CParty

#### Description:
Software implementation of CParty.     
CParty is an algorithm for computing the conditional partition function for density-2 RNA pseudoknots.

#### Cite: 
Trinity, Luke, Sebastian Will, Yann Ponty, Ulrike Stege, and Hosna Jabbari. "CParty: Conditional partition function for density-2 RNA pseudoknots." bioRxiv (2023): 2023-05.
https://doi.org/10.1101/2023.05.16.541023 

#### Supported OS: 
Linux 
macOS 

### Installation:  
Requirements: A compiler that supports C++11 standard (tested with g++ version 4.7.2 or higher)  and CMake version 3.1 or greater.    

[CMake](https://cmake.org/install/) version 3.1 or greater must be installed in a way that HFold can find it.    
To test if your Mac or Linux system already has CMake, you can type into a terminal:      
```
cmake --version
```
If it does not print a cmake version greater than or equal to 3.1, you will have to install CMake depending on your operating system.

#### Mac:    
Easiest way is to install homebrew and use that to install CMake.    
To do so, run the following from a terminal to install homebrew:      
```  
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"   
```    
When that finishes, run the following from a terminal to install CMake.     
```   
brew install cmake   
``` 
#### Linux:    
Run from a terminal     
```
wget http://www.cmake.org/files/v3.8/cmake-3.8.2.tar.gz
tar xzf cmake-3.8.2.tar.gz
cd cmake-3.8.2
./configure
make
make install
```
[Linux instructions source](https://geeksww.com/tutorials/operating_systems/linux/installation/downloading_compiling_and_installing_cmake_on_linux.php)

#### Steps for installation   
1. [Download the repository](https://github.com/HosnaJabbari/HFold.git) and extract the files onto your system.
2. From a command line in the root directory (where this README.md is) run
```
cmake -H. -Bbuild
cmake --build build
```   
If you need to specify a specific compiler, such as g++, you can instead run something like   
```
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
```   
This can be useful if you are getting errors about your compiler not having C++17 features.

Help
========================================

```
Usage: CParty[options] [input sequence]
```

Read input file from cmdline; predict minimum free energy, ensemble energy, and optimum structure using the RNA folding algorithm.


```
  -h, --help             Print help and exit
  -V, --version          Print version and exit
  -r, --input-structure  Give a restricted structure as an input structure
  -i, --input-file       Give a path to an input file containing the sequence (and input structure if known)
  -o, --output-file      Give a path to an output file which will the sequence, and its structure and energy
  -n, --opt              Specify the number of suboptimal structures to output (default is 1)
  -p  --pk-free          Specify whether you only want the pseudoknot-free structure to be calculated
  -k  --pk-only          Only add base pairs which cross the constraint structure. The constraint structure is returned if there are no energetically favorable crossing base pairs
  -d  --dangles          Specify the dangle model to be used (base is 2)
  -P, --paramFile        Read energy parameters from paramfile, instead of using the default parameter set.\n
      --noConv           Do not convert DNA into RNA. This will use the Matthews 2004 parameters for DNA
  
```


#### How to use:

    Remarks:
        make sure the <arguments> are enclosed in "", for example -r "..(...).." instead of -r ..(...)..
        input file for -i must be .txt
        if -i is provided with just a file name without a path, it is assuming the file is in the diretory where the executable is called
        if -o is provided with just a file name without a path, the output file will be generated in the diretory where the executable is called
        if -o is provided with just a file name without a path, and if -i is provided, then the output file will be generated in the directory where the input file is located
        if suboptimal structures are specified, repeated structures are skipped. That is, if different input structures come to the same conclusion, only those that are different are shown
        If no input structure is given, or suboptimal structures are greater than the number given, CParty generates hotspots to be used as input structures -- where hotspots are energetically favorable stems
        The default parameter file is Turner2004. This can be changed via -P and specifying the parameter file you would like
    
    Sequence requirements:
        containing only characters GCAU

    Structure requirements:
        -pseudoknot free
        -containing only characters .x()
        Remarks:
            Restricted structure symbols:
                () restricted base pair
                . no restriction
                x restricted to unpaired


    Input file requirements:
            Line1: FASTA name (optional)
            Line2: Sequence
            Line3: Structure
        sample:
            >Sequence1 (optional)
            GCAACGAUGACAUACAUCGCUAGUCGACGC
            (............................)

#### Example:
    Assume you are in the CParty directory
    ./build/CParty -i "/home/username/Desktop/myinputfile.txt"
    ./build/CParty -i "/home/username/Desktop/myinputfile.txt" --o "outputfile.txt"
    ./build/CParty -r "(............................)" GCAACGAUGACAUACAUCGCUAGUCGACGC
    ./build/CParty -r "(((((.........................)))))................" -d2 GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC
    ./build/CParty -p -r "(............................)" -o "/home/username/Desktop/some_folder/outputfile.txt" GCAACGAUGACAUACAUCGCUAGUCGACGC
    ./build/CParty -n 3 -r "(............................)" -o "/home/username/Desktop/some_folder/outputfile.txt" GCAACGAUGACAUACAUCGCUAGUCGACGC
    ./build/CParty -k -r "(............................)" GCAACGAUGACAUACAUCGCUAGUCGACGC
    ./build/CParty -P "src/params/parameters_DP09.txt" -r "(............................)" GCAACGAUGACAUACAUCGCUAGUCGACGC



### SARS-CoV-2 Example
    ./build/CParty -r "..(((((((((((..........)))))))))))....................................." UUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCU

    
## Questions
For questions, you can email mateo2@ualberta.ca