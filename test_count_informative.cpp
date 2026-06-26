#include "mp_features_export.h"
#include "alignment.h"
#include "tools.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <alignment_file>" << std::endl;
        return 1;
    }
    
    // Load alignment
    Params params;
    params.aln_file = argv[1];
    
    Alignment* aln = new Alignment(params.aln_file, params.sequence_type, params.intype);
    if (!aln) {
        std::cerr << "ERROR: Failed to load alignment from " << params.aln_file << std::endl;
        return 1;
    }
    
    // Count informative sites
    int informative_count = countInformativeSites(aln);
    
    // Print results
    std::cout << "Alignment: " << params.aln_file << std::endl;
    std::cout << "Total sites: " << aln->getNSite() << std::endl;
    std::cout << "Total patterns: " << aln->getNPattern() << std::endl;
    std::cout << "Number of taxa: " << aln->getNSeq() << std::endl;
    std::cout << "Number of states: " << aln->num_states << std::endl;
    std::cout << "Informative sites: " << informative_count << std::endl;
    std::cout << "Non-informative sites: " << (aln->getNSite() - informative_count) << std::endl;
    
    delete aln;
    return 0;
}

// Dummy implementations to satisfy linker for functions defined in pda.cpp
void printCopyrightMP(std::ostream &out) {}
void summarizeHeader(std::ostream &out, Params &params, bool b, InputType intype) {}
void summarizeFooter(std::ostream &out, Params &params) {}

class PDNetwork;
class SplitSet;
struct PDRelatedMeasures;

void summarizeSplit(Params& params, PDNetwork& net, std::vector<SplitSet>& ss, PDRelatedMeasures& pdrm, bool b) {}
void printCopyright(std::ostream& out) {}

