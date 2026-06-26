/*
 * mp_features_export.h
 * 
 * MP Features Export Module for MPBoot
 * 
 * This module provides functionality to export parsimony computation data
 * for machine learning-based phylogenetic difficulty prediction.
 * 
 * Usage:
 *   mpboot --export-mp-features -s alignment.fasta -z trees.nwk -pre output
 * 
 * This generates output.mp_features.json containing:
 *   - Tree scores for K=20 MP trees
 *   - Pattern-level parsimony scores
 *   - Pattern frequencies
 *   - Alignment metadata
 */

#ifndef MP_FEATURES_EXPORT_H
#define MP_FEATURES_EXPORT_H

#include <vector>
#include <string>

// Forward declarations
class Alignment;
class ParsTree;
struct Params;

/**
 * Data structure holding parsimony feature export data
 * This is the C++ representation that gets serialized to JSON
 */
struct MPFeaturesData {
    // K=20 tree scores (one integer per tree)
    std::vector<int> tree_scores;
    
    // Pattern-level scores: [K trees × N patterns]
    // Stored as flattened 2D array: pattern_scores[tree_idx * num_patterns + pattern_idx]
    std::vector<int> pattern_scores;
    
    // Pattern frequencies (N values, one per pattern)
    // These map patterns back to sites for site-wise statistics
    std::vector<int> pattern_freqs;
    
    // Alignment metadata
    int num_sites;              // Total number of sites in alignment (getNSite())
    int num_patterns;           // Number of unique patterns (size())
    int num_informative_sites;  // Sites meeting parsimony informativeness criteria
    int num_taxa;               // Number of sequences (getNSeq())
    int num_states;             // Number of character states (num_states: 4=DNA, 20=protein)
    
    MPFeaturesData() : num_sites(0), num_patterns(0), 
                       num_informative_sites(0), num_taxa(0), num_states(0) {}
};

/**
 * Main entry point for MP features export mode
 * 
 * @param params Global parameters structure
 * @param aln Alignment pointer
 * @param tree_strings Vector of Newick tree strings from candidate trees
 * @return 0 on success, non-zero error code on failure
 */
int exportMPFeatures(Params &params, Alignment *aln, const std::vector<std::string> &tree_strings);

/**
 * Load K trees from treefile into vector
 * 
 * @param treefile_path Path to Newick format treefile
 * @param aln Alignment pointer (for taxon name validation)
 * @param trees (OUT) Vector of ParsTree pointers
 * @param max_trees Maximum number of trees to read (default -1 means read all)
 * @return true on success, false on error
 */
bool loadTrees(const char *treefile_path, Alignment *aln, std::vector<ParsTree*> &trees, int max_trees = -1);

/**
 * Compute parsimony data for all K trees
 * 
 * @param trees Vector of K ParsTree objects
 * @param aln Alignment pointer
 * @param data (OUT) Filled MPFeaturesData structure
 */
void computeParsimonyData(std::vector<ParsTree*> &trees, Alignment *aln, MPFeaturesData &data);

/**
 * Count parsimony informative sites in alignment
 * A site is informative if it has at least 2 different states, 
 * each appearing in at least 2 taxa
 * 
 * @param aln Alignment pointer
 * @return count of informative sites
 */
int countInformativeSites(Alignment *aln);

/**
 * Serialize MPFeaturesData to JSON file
 * 
 * @param data Parsimony data structure
 * @param output_path Path to output JSON file
 * @return true on success, false on I/O error
 */
bool serializeToJSON(const MPFeaturesData &data, const char *output_path);

#endif // MP_FEATURES_EXPORT_H
