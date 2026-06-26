#include "parstree.h"
#include "mp_features_export.h"
#include "alignment.h"
#include "tools.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <alignment_file> <treefile>" << std::endl;
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
    
    std::cout << "Alignment loaded: " << aln->getNSeq() << " taxa, " 
              << aln->getNSite() << " sites" << std::endl;
    
    // Load trees
    std::vector<ParsTree*> trees;
    bool success = loadTrees(argv[2], aln, trees);
    
    if (!success) {
        std::cerr << "ERROR: Failed to load trees" << std::endl;
        delete aln;
        return 1;
    }
    
    std::cout << "\nSuccessfully loaded " << trees.size() << " tree(s)" << std::endl;
    
    // Print tree info
    for (size_t i = 0; i < trees.size(); i++) {
        std::cout << "Tree " << (i+1) << ": " 
                  << trees[i]->leafNum << " leaves, "
                  << trees[i]->nodeNum << " total nodes" << std::endl;
    }
    
    // Clean up
    for (auto t : trees) delete t;
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
