/*
 * mp_features_export.cpp
 * 
 * Implementation of MP Features Export Module for MPBoot
 */

#include "mp_features_export.h"
#include "phylotree.h"
#include "alignment.h"
#include "parstree.h"
#include "tools.h"
#include "timeutil.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

bool parseTreeStrings(const vector<string> &tree_strings, Alignment *aln, vector<ParsTree*> &trees) {
    trees.clear();
    for (size_t i = 0; i < tree_strings.size(); i++) {
        try {
            ParsTree *tree = new ParsTree(aln);
            tree->readTreeString(tree_strings[i]);
            
            // Basic taxon count validation
            if (tree->leafNum != aln->getNSeq()) {
                cerr << "ERROR: Parsed tree " << (i+1) << " has " << tree->leafNum 
                     << " taxa, but alignment has " << aln->getNSeq() << endl;
                delete tree;
                for (auto t : trees) delete t;
                trees.clear();
                return false;
            }
            trees.push_back(tree);
        } catch (const char* str) {
            cerr << "ERROR: Failed to parse tree string " << (i+1) << ": " << str << endl;
            for (auto t : trees) delete t;
            trees.clear();
            return false;
        }
    }
    return true;
}

/**
 * Main entry point for MP features export mode
 */
int exportMPFeatures(Params &params, Alignment *aln, const vector<string> &tree_strings) {
    if (!params.catphy_mode) {
        cerr << "ERROR: --export-mp-features requires CatPhy mode to be enabled (-catphy)" << endl;
        return 1;
    }
    
    if (!params.out_prefix) {
        cerr << "ERROR: --export-mp-features requires -pre (output prefix)" << endl;
        return 1;
    }

    cout << "Starting MP Features Export Module using " << tree_strings.size() << " candidate trees..." << endl;
    
    // Parse the in-memory tree strings
    vector<ParsTree*> trees;
    if (!parseTreeStrings(tree_strings, aln, trees)) {
        cerr << "ERROR: Failed to parse candidate tree strings for feature extraction." << endl;
        return 1;
    }
    
    cout << "Successfully parsed " << trees.size() << " candidate trees for parsimony computation." << endl;
    
    // Create data structure and compute parsimony
    MPFeaturesData data;
    double start_time = getRealTime();
    computeParsimonyData(trees, aln, data);
    cout << "Parsimony computation completed in " << (getRealTime() - start_time) << " seconds." << endl;
    
    // Write out the candidate trees to a .treefile so Python can compute topology features
    string treefile_path = string(params.out_prefix) + ".mp_features.treefile";
    ofstream tree_out(treefile_path.c_str());
    if (tree_out.is_open()) {
        // Add a standard Newick comment tag at the top of the file so users 
        // quickly know how many trees are inside without counting lines.
        tree_out << "[ " << tree_strings.size() << " candidate parsimony trees exported by MPBoot ]" << endl;
        
        for (const string &str : tree_strings) {
            tree_out << str << endl;
        }
        tree_out.close();
    } else {
        cerr << "ERROR: Failed to write trees to " << treefile_path << endl;
    }
    
    // Determine JSON output file path
    string output_path = string(params.out_prefix) + ".mp_features.json";
    
    // Serialize to JSON
    if (!serializeToJSON(data, output_path.c_str())) {
        cerr << "ERROR: Failed to serialize features to JSON." << endl;
        // Clean up allocated trees before returning
        for (auto t : trees) delete t;
        return 1;
    }
    
    // Clean up
    for (auto t : trees) delete t;
    
    cout << "Successfully exported MP features to " << output_path << endl;
    return 0;
}

/**
 * Load K trees from treefile into vector
 */
bool loadTrees(const char *treefile_path, Alignment *aln, vector<ParsTree*> &trees, int max_trees) {
    // Open treefile as input stream
    ifstream in(treefile_path);
    if (!in.is_open()) {
        cerr << "ERROR: Cannot open treefile: " << treefile_path << endl;
        return false;
    }
    
    cout << "Reading trees from " << treefile_path << " ..." << endl;
    if (max_trees > 0) {
        cout << "  (Will read up to " << max_trees << " trees)" << endl;
    }
    
    // Clear any existing trees
    trees.clear();
    
    // Loop to read trees until EOF or error
    int tree_count = 0;
    while (true) {
        if (max_trees > 0 && tree_count >= max_trees) {
            break;
        }
        // Skip whitespace manually to avoid MTree::readTree throwing EOF errors
        // because it calls outError() and exit(1) internally instead of returning
        char ch;
        while (in.get(ch)) {
            if (ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n') {
                in.putback(ch);
                break;
            }
        }
        
        if (in.eof()) {
            break; // Valid EOF reached
        }

        try {
            // Create new ParsTree object
            ParsTree *tree = new ParsTree();
            
            // Read tree using MTree::readTree()
            bool is_rooted = false;
            tree->readTree(in, is_rooted);
            
            // Check if tree was successfully read
            if (tree->leafNum == 0) {
                // Empty tree, probably EOF
                delete tree;
                break;
            }
            
            // Validate taxon names match alignment
            if (tree->leafNum != aln->getNSeq()) {
                cerr << "ERROR: Tree " << (tree_count + 1) 
                     << " has " << tree->leafNum << " taxa, but alignment has " 
                     << aln->getNSeq() << " taxa" << endl;
                delete tree;
                // Clean up already loaded trees
                for (auto t : trees) delete t;
                trees.clear();
                in.close();
                return false;
            }
            
            // Validate taxon names
            bool names_match = true;
            for (int i = 0; i < tree->leafNum; i++) {
                string tree_taxon = tree->findNodeID(i)->name;
                string aln_taxon = aln->getSeqName(i);
                
                // Check if this taxon exists in alignment
                bool found = false;
                for (int j = 0; j < aln->getNSeq(); j++) {
                    if (tree_taxon == aln->getSeqName(j)) {
                        found = true;
                        break;
                    }
                }
                
                if (!found) {
                    cerr << "ERROR: Tree " << (tree_count + 1) 
                         << " contains taxon '" << tree_taxon 
                         << "' not found in alignment" << endl;
                    names_match = false;
                    break;
                }
            }
            
            if (!names_match) {
                delete tree;
                // Clean up already loaded trees
                for (auto t : trees) delete t;
                trees.clear();
                in.close();
                return false;
            }
            
            // Add tree to vector
            trees.push_back(tree);
            tree_count++;
            
        } catch (const char* str) {
            cerr << "ERROR: Failed to read tree " << (tree_count + 1) 
                 << " from " << treefile_path << ": " << str << endl;
            // Clean up already loaded trees
            for (auto t : trees) delete t;
            trees.clear();
            in.close();
            return false;
        } catch (...) {
            // Might be EOF or other error
            break;
        }
    }
    
    in.close();
    
    cout << "Successfully loaded " << tree_count << " tree(s)" << endl;
    
    if (tree_count == 0) {
        cerr << "ERROR: No trees found in " << treefile_path << endl;
        return false;
    }
    
    return true;
}

/**
 * Compute parsimony data for all K trees
 */
void computeParsimonyData(vector<ParsTree*> &trees, Alignment *aln, MPFeaturesData &data) {
    int K = trees.size();
    int N = aln->getNPattern();
    
    // Initialize metadata from alignment
    data.num_patterns = N;
    data.num_sites = aln->getNSite();
    data.num_taxa = aln->getNSeq();
    data.num_states = aln->num_states;
    data.num_informative_sites = countInformativeSites(aln);
    
    // Allocate arrays
    data.tree_scores.resize(K);
    data.pattern_scores.resize(K * N);  // Flattened 2D array
    data.pattern_freqs.resize(N);
    
    // Extract pattern frequencies (only need to do this once)
    for (int ptn = 0; ptn < N; ptn++) {
        data.pattern_freqs[ptn] = aln->at(ptn).frequency;
    }
    
    cout << "Computing parsimony scores for " << K << " trees..." << endl;
    
    // For each tree: compute parsimony scores
    for (int t = 0; t < K; t++) {
        ParsTree *tree = trees[t];
        
        // Set alignment
        tree->setAlignment(aln);
        
        // Initialize cost matrix for Fitch parsimony (uniform cost)
        // This must be done before computeParsimony()
        tree->loadCostMatrixFile((char*)"fitch");
        
        // Set root node (use first taxon as root)
        tree->setRootNode((char*)aln->getSeqName(0).c_str());
        
        // Fix degree-2 nodes (often created by rooted trees in input) which break 
        // the strictly-bifurcating assumption in ParsTree::computePartialParsimony
        NodeVector internal_nodes;
        tree->getInternalNodes(internal_nodes);
        for (int i = 0; i < internal_nodes.size(); i++) {
            PhyloNode* node = (PhyloNode*)internal_nodes[i];
            if (node && node->degree() == 2 && node != tree->root) {
                Neighbor* nei1 = node->neighbors[0];
                Neighbor* nei2 = node->neighbors[1];
                Node* child1 = nei1->node;
                Node* child2 = nei2->node;
                double new_len = nei1->length + nei2->length;
                
                child1->updateNeighbor(node, child2, new_len);
                child2->updateNeighbor(node, child1, new_len);
                
                node->neighbors.clear();
            }
        }
        
        // Initialize partial parsimony vectors
        tree->initializeAllPartialPars();
        
        // Compute parsimony score for this tree
        int tree_score = tree->computeParsimony();
        data.tree_scores[t] = tree_score;
        
        // Copy pattern scores from _pattern_pars array
        // _pattern_pars is filled by computeParsimony()
        for (int ptn = 0; ptn < N; ptn++) {
            data.pattern_scores[t * N + ptn] = tree->getPatternPars(ptn);
        }
        
        if ((t + 1) % 5 == 0 || t == K - 1) {
            cout << "  Processed " << (t + 1) << "/" << K << " trees" << endl;
        }
    }
    
    cout << "Parsimony computation complete." << endl;
    cout << "  Tree scores range: [" << *min_element(data.tree_scores.begin(), data.tree_scores.end())
         << ", " << *max_element(data.tree_scores.begin(), data.tree_scores.end()) << "]" << endl;
}

/**
 * Count parsimony informative sites in alignment
 * A site is informative if it has at least 2 different states,
 * each appearing in at least 2 taxa
 */
int countInformativeSites(Alignment *aln) {
    int informative_count = 0;
    int num_patterns = aln->getNPattern();
    int num_taxa = aln->getNSeq();
    int num_states = aln->num_states;
    char state_unknown = aln->STATE_UNKNOWN;
    
    // Iterate through each pattern
    for (int ptn = 0; ptn < num_patterns; ptn++) {
        Pattern &pattern = aln->at(ptn);
        
        // Count frequency of each state in this pattern
        // Only count standard (unambiguous) states: 0 to num_states-1
        int state_freq[256] = {0};
        for (int seq = 0; seq < num_taxa; seq++) {
            char state = pattern[seq];
            
            // Only count unambiguous states (not gaps, not unknown, not ambiguous)
            // Standard states are in range [0, num_states-1]
            if (state >= 0 && state < num_states) {
                state_freq[(unsigned char)state]++;
            }
        }
        
        // Count how many states appear at least 2 times
        int states_with_min_2 = 0;
        for (int i = 0; i < num_states; i++) {
            if (state_freq[i] >= 2) {
                states_with_min_2++;
            }
        }
        
        // Site is informative if at least 2 states appear >= 2 times
        if (states_with_min_2 >= 2) {
            // Add frequency of this pattern to count (patterns can represent multiple sites)
            informative_count += pattern.frequency;
        }
    }
    
    return informative_count;
}

/**
 * Serialize MPFeaturesData to JSON file
 */
bool serializeToJSON(const MPFeaturesData &data, const char *output_path) {
    ofstream out(output_path);
    if (!out.is_open()) {
        cerr << "ERROR: Cannot write to file: " << output_path << endl;
        return false;
    }
    
    // Write JSON with 2-space indentation
    out << "{\n";
    
    // Tree scores array
    out << "  \"tree_scores\": [";
    for (size_t i = 0; i < data.tree_scores.size(); i++) {
        if (i > 0) out << ", ";
        out << data.tree_scores[i];
    }
    out << "],\n";
    
    // Pattern scores array (flattened 2D)
    out << "  \"pattern_scores\": [";
    for (size_t i = 0; i < data.pattern_scores.size(); i++) {
        if (i > 0) out << ", ";
        out << data.pattern_scores[i];
    }
    out << "],\n";
    
    // Pattern frequencies array
    out << "  \"pattern_freqs\": [";
    for (size_t i = 0; i < data.pattern_freqs.size(); i++) {
        if (i > 0) out << ", ";
        out << data.pattern_freqs[i];
    }
    out << "],\n";
    
    // Metadata fields
    out << "  \"num_sites\": " << data.num_sites << ",\n";
    out << "  \"num_patterns\": " << data.num_patterns << ",\n";
    out << "  \"num_informative_sites\": " << data.num_informative_sites << ",\n";
    out << "  \"num_taxa\": " << data.num_taxa << ",\n";
    out << "  \"num_states\": " << data.num_states << "\n";
    
    out << "}\n";
    
    out.close();
    
    cout << "JSON output written to " << output_path << endl;
    return true;
}
