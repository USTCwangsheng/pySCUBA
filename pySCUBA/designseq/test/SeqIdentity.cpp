/*
 * SeqIdentity.cpp
 *
 *  Created on: 2019/3/21
 *      Author: hbfrank
 */

#include <iostream>
#include <string>

void calculateSeqIdentity(const std::string & first, const std::string & second)
{
    if (first.size() != second.size()) {
        std::cerr << "[Error] The sequences are not the same length." << std::endl;
        exit(1);
    }
    if (first.empty()) {
        std::cerr << "[Error] Empty sequences." << std::endl;
        exit(1);
    }
    int numSites = first.size();
    int numIdenticalSites = 0;
    for (int i = 0; i < numSites; i++) {
        if (first[i] == second[i]) {
            numIdenticalSites++;
        }
    }
    double identity = ((double)numIdenticalSites) / numSites;
    std::cout << numSites << " " << numIdenticalSites << " " << identity << std::endl;
}

void printUsage(const std::string & selfname)
{
    std::cout << "Usage:" << std::endl;
    std::cout << "$ " << selfname << " <firstSequence> <secondSequence>"
            << std::endl;
}

int main(int argc, char** argv)
{
    std::string selfname = argv[0];

    if (argc != 3) {
        std::cerr << "[Error] Missing arguments." << std::endl;
        printUsage(selfname);
        return 1;
    }

    std::string firstSeq(argv[1]);
    std::string secondSeq(argv[2]);

    calculateSeqIdentity(firstSeq, secondSeq);

    return 0;
}

