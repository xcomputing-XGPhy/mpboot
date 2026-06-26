#include "parstree.h"
#include "mp_features_export.h"
#include "alignment.h"
#include "tools.h"
#include <iostream>
#include <vector>

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
        std::cerr << "ERROR: Failed to load alignment" << std::endl;
        return 1;
    }
    
    std::cout << "Alignment: " << aln->getNSeq() << " taxa, " 
              << aln->getNSite() << " sites, "
              << aln->getNPattern() << " patterns" << std::endl;
    
    // Load trees
    std::vector<ParsTree*> trees;
    if (!loadTrees(argv[2], aln, trees)) {
        std::cerr << "ERROR: Failed to load trees" << std::endl;
        delete aln;
        return 1;
    }
    
    std::cout << "Loaded " << trees.size() << " trees\n" << std::endl;
    
    // Compute parsimony data
    MPFeaturesData data;
    computeParsimonyData(trees, aln, data);
    
    // Print results
    std::cout << "\n=== Parsimony Data Summary ===" << std::endl;
    std::cout << "Metadata:" << std::endl;
    std::cout << "  num_sites: " << data.num_sites << std::endl;
    std::cout << "  num_patterns: " << data.num_patterns << std::endl;
    std::cout << "  num_informative_sites: " << data.num_informative_sites << std::endl;
    std::cout << "  num_taxa: " << data.num_taxa << std::endl;
    std::cout << "  num_states: " << data.num_states << std::endl;
    
    std::cout << "\nTree Scores:" << std::endl;
    for (size_t i = 0; i < data.tree_scores.size(); i++) {
        std::cout << "  Tree " << (i+1) << ": " << data.tree_scores[i] << std::endl;
    }
    
    std::cout << "\nPattern Frequencies (first 10):" << std::endl;
    for (size_t i = 0; i < std::min(data.pattern_freqs.size(), size_t(10)); i++) {
        std::cout << "  Pattern " << i << ": " << data.pattern_freqs[i] << std::endl;
    }
    
    std::cout << "\nPattern Scores (first tree, first 10 patterns):" << std::endl;
    for (size_t i = 0; i < std::min(data.num_patterns, 10); i++) {
        std::cout << "  Pattern " << i << ": " << data.pattern_scores[i] << std::endl;
    }
    
    // Verify dimensions
    std::cout << "\n=== Data Validation ===" << std::endl;
    std::cout << "tree_scores size: " << data.tree_scores.size() 
              << " (expected: " << trees.size() << ")" << std::endl;
    std::cout << "pattern_scores size: " << data.pattern_scores.size() 
              << " (expected: " << (trees.size() * data.num_patterns) << ")" << std::endl;
    std::cout << "pattern_freqs size: " << data.pattern_freqs.size() 
              << " (expected: " << data.num_patterns << ")" << std::endl;
    
    // Clean up
    for (auto t : trees) delete t;
    delete aln;
    
    std::cout << "\nTest completed successfully!" << std::endl;
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
