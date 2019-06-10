//
//  main.cpp
//  greedRelate
//  Greedy approach of keeping samples
//  Created by Shing Wan Choi on 04/06/2018.
//  Copyright © 2018 Shing Wan Choi. All rights reserved.
//

#include "misc.hpp"
#include <algorithm>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Useful for debugging
bool is_problem(const std::string& input)
{
    std::vector<std::string> prob = {};
    for (auto&& p : prob) {
        if (input == p) return true;
    }
    return false;
}

class Sample
{
public:
    // First occurance = itself, therefore occurance = 0
    Sample(std::string name, double pheno, double rand)
        : m_name(name), m_phenotype(pheno), m_rand_number(rand)
    {
        m_occur = 0;
    }
    // add a relative to the current sample
    int add(Sample* related)
    {
        // not sure why
        // TODO: See if we can skip the step of pushing nullptr
        // m_relatives.push_back(nullptr);
        // m_relatives.back() = related;
        m_relatives.push_back(related);
        m_removed = false;
        m_occur++;
        return m_occur;
    }
    int remove(std::ostream& os)
    {
        // First, check if we should remove any relatives before ourselves
        // If there are any to be removed, remove them first.Then
        // check if we still need to remove ourselves
        // If we no longer need to remove ourselves, immediately stop
        // transversing the relative vector
        // Do sorting first to ensure we remove the most related pair first
        std::sort(m_relatives.begin(), m_relatives.end(),
                  Sample::compare_sample);
        for (auto&& relative : m_relatives) {
            if (relative->removed() || relative->m_occur < m_occur
                || relative->m_phenotype > m_phenotype
                || relative->m_rand_number < m_rand_number
                || relative->m_name < m_name)
            {
                // do nothing
            }
            else if (!relative->removed())
            {

                // relative->m_occur > m_occur
                // relative->m_occur == m_occur
                //      relative->m_phenotype < m_phenotype
                //      relative->m_phenotype == m_phenotype
                //          relative->m_rand_number < m_rand_number
                //          relative->m_rand_number == m_rand_number
                //              relative->m_name < m_name
                // rescued by relative
                // try to remove the relative now
                relative->remove(os);
                // If we have nothing left, end.
                // If we can still remove, then continue on
                if (m_occur <= 0 || m_removed) return 0;
            }
        }
        // if we reach here, it means this sample still need to be removed.
        os << m_name << "\t"<< m_name <<  "\t" << m_occur << std::endl;
        m_occur = -1;
        m_removed = true;
        // Update our relatives to inform them we are dead
        for (auto&& relative : m_relatives) {
            relative->m_occur--;
            if (relative->m_occur <= 0) relative->m_removed = true;
        }
        // Now update our relative's ordering before removal. Most likely
        // not required. But should be safer this way as the relative vector
        // should be tiny (< 100 ) anyway.
        std::sort(m_relatives.begin(), m_relatives.end(),
                  Sample::compare_sample);
        for (auto&& relative : m_relatives) {
            if (!relative->removed() && relative->m_occur > 0) {
                relative->remove(os);
            }
        }
        return m_occur;
    }

    std::string debug() const
    {
        std::string occur;          // string which will contain the result
        std::ostringstream convert; // stream used for the conversion
        convert << m_occur; // insert the textual representation of 'Number' in
                            // the characters in the stream
        occur = convert.str();
        return m_name + " " + occur;
    }

    static bool compare_sample(Sample const* a, Sample const* b)
    {
        // Ordering is based on
        // 1. Who has more related pair? More is first
        // 2. Who has a larger phenotype? Higher come later (reserve cases)
        // 3. Who has a higher random number? Higher come first
        // 4. Who has a larger ID? Do string comparison, larger come first
        if (a->m_occur == b->m_occur) {
            if (misc::logically_equal(a->m_phenotype, b->m_phenotype)) {
                if (misc::logically_equal(a->m_rand_number, b->m_rand_number)) {
                    if (a->m_name == b->m_name)
                        return false;
                    else
                        return a->m_name > b->m_name;
                }
                else
                    return a->m_rand_number > b->m_rand_number;
            }
            else
                return a->m_phenotype < b->m_phenotype;
        }
        else
            return a->m_occur > b->m_occur;
    }

    int get_occur() const { return m_occur; }
    std::string get_name() const { return m_name; }
    double get_pheno() const { return m_phenotype; }
    double get_rand() const { return m_rand_number; }
    bool removed() const
    {
        // add additional condition where the current sample is removed
        return (m_occur == 0 || m_removed);
    }

private:
    std::vector<Sample*> m_relatives;
    std::string m_name;
    double m_phenotype;
    double m_rand_number; // This is a random number to solve the tie
    int m_occur;
    bool m_removed = true; // so any samples with invalid pairs will be ignoreds
};

void usage()
{
    fprintf(stderr, " GreedyRelate\n");
    fprintf(stderr, " Sam Choi\n");
    fprintf(stderr, " v1.2.0 ( 2019-06-04 )\n");
    fprintf(stderr, " ==============================\n");
    fprintf(stderr, " This programme will try to minize the number of samples "
                    "that\n");
    fprintf(stderr, "need to be removed due to relatedness\n");
    fprintf(stderr,
            " Samples that should be removed are output to the STDOUT or\n ");
    fprintf(stderr, "the Output file specified by --out");
    fprintf(stderr, " \n");
    fprintf(stderr, " Usage: GreedyRelate [options] -r <relatedness>\n");
    fprintf(stderr, "       -r | --relate      Relationship file (Required)\n");
    fprintf(stderr, "       -p | --pheno       Phenotype file\n");
    fprintf(stderr, "       -t | --threshold   Relatedness Threshold\n");
    fprintf(stderr, "       -k | --keep        Ignore any samples not presented\n");
    fprintf(stderr, "                          in this file. We only use the first\n");
    fprintf(stderr, "                          column as the ID.\n");
    fprintf(stderr, "       -i | --id1         Column containing the first ID\n");
    fprintf(stderr, "                          When provided, will assume PLINK\n");
    fprintf(stderr, "                          like formats\n");
    fprintf(stderr, "       -I | --id2         Column containing the second ID\n");
    fprintf(stderr, "                          When provided, will assume PLINK\n");
    fprintf(stderr, "                          like formats\n");
    fprintf(stderr, "       -f | --fstat       Column containing the F stat\n");
    fprintf(stderr, "                          When provided, will assume PLINK\n");
    fprintf(stderr, "                          like formats\n");
    fprintf(stderr, "       -P | --plink       Input is a plink format. \n");
    fprintf(stderr, "                          Will set -i IID1 -I IID2 -f PI_HAT\n");
    fprintf(stderr, "                          Will over-ride -i -I and -f\n");
    fprintf(stderr,
            "       -o | --out         Output name. Stdout if not provided\n");
    fprintf(stderr,
            "       -s | --seed        Seed for the random number generator\n");
    fprintf(stderr,
            "       -h | --help        Display this help message\n\n\n");
    fprintf(stderr, " Details:\n");
    fprintf(stderr,
            "       Relationship file should have the following format:\n");
    fprintf(stderr, "       ID    Pair    Factor\n\n");
    fprintf(stderr, "       We do assume there is a header line\n\n");
    fprintf(stderr,
            "       Phenotype file should have the following format:\n");
    fprintf(stderr, "       ID    Pheno\n\n");
    fprintf(stderr, "       Again, we do assume there is a header line\n");
    fprintf(stderr,
            "       Note: Phenotype information only used to decide which\n");
    fprintf(stderr,
            "             samples to leave behind when there is a tie,\n");
    fprintf(stderr,
            "             where sample with higher phenotype value will be\n");
    fprintf(stderr, "             retained.\n");
    fprintf(stderr,
            "             When no phenotype information is provided, \n");
    fprintf(stderr,
            "             we will randomly select one sample to remove\n");
}


template <typename T>
void delete_pointed_to(T* const ptr)
{
    delete ptr;
}


std::vector<Sample*> kin3col(const std::string &relate_name,
                             const std::unordered_set<std::string> &include_samples,
                             const std::unordered_map<std::string, double> &phenotype,
                             const bool &keep_samples, const double &threshold,
                             const std::mt19937 &rand_gen){
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::unordered_map<std::string, size_t> sample_index;
    std::unordered_map<size_t, size_t> direction; // First size_t = pair, second
                                                  // size_t = index of the
                                                  // neighbour
    std::vector<Sample*> sample_list;
    std::ifstream relate;
    relate.open(relate_name.c_str());
    if (!relate.is_open()) {
        fprintf(stderr, "ERROR: Cannot open relationship file %s\n",
                relate_name.c_str());
        exit(-1);
    }

    std::string line;
    // Assume there is a header
    // Or we can add in
    fprintf(stderr,
            "Assuming there is a header file for the relatedness file\n");
    getline(relate, line);
    std::vector<std::string> token;
    size_t sample_idx = 0;
    std::unordered_set<size_t> remove_pair;
    while (getline(relate, line)) {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() != 3) {
            fprintf(stderr, "ERROR: Relationship file format incorrect! "
                            "Require 3 columns\n");
            exit(-1);
        }
        std::string id = token[0];
        size_t pair = 0;
        double factor = 0.0;
        try
        {
            pair = misc::convert<size_t>(token[1]);
            factor = misc::convert<double>(token[2]);
        }
        catch (const std::runtime_error&)
        {
            fprintf(stderr, "ERROR: Cannot convert some of the information in "
                            "the relationship file\n");
            fprintf(stderr, "Input: %s\n", line.c_str());
            exit(-1);
        }
        // Now we have id pair and factor
        // Ignore any pairs with relatedness less than our threshold
        if (factor <= threshold) continue;
        if(keep_samples && include_samples.find(id) == include_samples.end()){
            remove_pair.insert(pair);
            continue;
        }
        double pheno = -9;
        if (phenotype.find(id) != phenotype.end()) pheno = phenotype.at(id);
        if (sample_index.find(id) == sample_index.end()) {
            // this is a new sample
            sample_list.push_back(
                new Sample(id, pheno, distribution(rand_gen)));
            sample_index[id] = sample_idx;
            ++sample_idx;
        }
        // worst case scenario in UKBB -> 109 relatives
        // still better than resorting the whole vector
        if (direction.find(pair) != direction.end()) {
            // this is not a new pair
            size_t dir_id = direction[pair];
            size_t sam_id = sample_index[id];
            sample_list[dir_id]->add(sample_list[sam_id]);
            sample_list[sam_id]->add(sample_list[dir_id]);
        }
        else
        {
            // a new pair is located
            direction[pair] = sample_index[id];
        }
    }
    return sample_list;
}

std::vector<Sample*> plink_format_process(const std::string &relate_name,
                                          const std::unordered_set<std::string> &include_samples,
                                          const std::unordered_map<std::string, double> &phenotype,
                                          const bool &keep_samples, const double &threshold,
                                          const std::mt19937 &rand_gen,
                                          const std::string &id_1_col,
                                          const std::string &id_2_col,
                                          const std::string &f_col){
    std::vector<Sample*> sample_list;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::unordered_map<std::string, size_t> sample_index;
    std::ifstream relate;
    relate.open(relate_name.c_str());
    if (!relate.is_open()) {
        fprintf(stderr, "ERROR: Cannot open relationship file %s\n",
                relate_name.c_str());
        exit(-1);
    }
    std::string line;
    // there must be a header
    std::getline(relate, line);
    std::vector<std::string> token;
    misc::trim(line);
    if(line.empty()){
        throw std::runtime_error("Erorr: First line of the related file cannot be empty!");
    }
    token=misc::split(line);
    // now get the index
    size_t id1_idx=0, id2_idx=0, f_idx=0;
    for(size_t i =  0; i < token.size(); ++i){
        if(token[i]==id_1_col) id1_idx=i;
        else if(token[i] ==id_2_col) id2_idx=i;
        else if(token[i]==f_col) f_idx=i;
        if(id1_idx && id2_idx && f_idx) break; // found all
    }
    size_t max_idx = id1_idx;
    if(id2_idx > max_idx) max_idx = id2_idx;
    if(f_idx > max_idx) max_idx = f_idx;
    while(std::getline(relate, line)){
        misc::trim(line);
        if(line.empty()) continue;
        token = misc::split(line);
        if(token.size() <= max_idx){
            fprintf(stderr, "ERROR: Relationship file format incorrect! "
                            "Require at least %lu columns\n", max_idx);
            exit(-1);
        }
        std::string id1 = token[id1_idx];
        std::string id2 = token[id2_idx];
        double fstat = 0;
        try {
            fstat = misc::convert<double>(token[f_idx]);
        } catch(const std::runtime_error&)
        {
            fprintf(stderr, "ERROR: Cannot convert some of the information in "
                            "the relationship file\n");
            fprintf(stderr, "Input: %s\n", line.c_str());
            exit(-1);
        }
        if (fstat <= threshold) continue;
        if(keep_samples && (include_samples.find(id1) == include_samples.end()||
               include_samples.find(id2) == include_samples.end()) ){
            continue;
        }
        double pheno1 = -9, pheno2=-9;
        if (phenotype.find(id1) != phenotype.end()) pheno1 = phenotype.at(id1);
        if(phenotype.find(id2)!= phenotype.end()) pheno2 = phenotype.at(id2);
        // first, check if we need to build a new sample
        auto id1_loc = sample_index.find(id1);
        auto id2_loc = sample_index.find(id2);
        if(id1_loc == sample_index.end()){
            sample_list.emplace_back(new Sample(id1, pheno1, distribution(rand_gen)));
            sample_index[id1] = sample_list.size()-1;
            id1_loc = sample_index.find(id1);
        }
        if(id2_loc == sample_index.end()){
            sample_list.emplace_back(new Sample(id2, pheno2, distribution(rand_gen)));
            sample_index[id2] = sample_list.size()-1;
            id1_loc = sample_index.find(id2);
        }
        sample_list[id1_loc->second]->add(sample_list[id2_loc->second]);
        sample_list[id2_loc->second]->add(sample_list[id1_loc->second]);
    }
    return sample_list;
}

int main(int argc, char* argv[])
{
    if (argc <= 1) {
        usage();
        exit(0);
    }
    int use_plink_default = false;
    static const char* optString = "r:p:t:s:n:o:k:i:I:f:h?";
    static const struct option longOpts[] = {
        {"relate", required_argument, nullptr, 'r'},
        {"pheno", required_argument, nullptr, 'p'},
        {"threshold", required_argument, nullptr, 't'},
        {"thread", required_argument, nullptr, 'n'},
        {"seed", required_argument, nullptr, 't'},
        {"out", required_argument, nullptr, 'o'},
    {"keep", required_argument, nullptr, 'k'},
    {"id1", required_argument, nullptr, 'i'},
    {"id2", required_argument, nullptr, 'I'},
    {"help", no_argument, &use_plink_default, 1},
    {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};
    std::string relate_name = "";
    std::string pheno_name = "";
    std::string out_name = "";
    std::string keep_name = "";
    std::string id_1_col = "";
    std::string id_2_col = "";
    std::string f_col = "";
    bool provide_seed = false;
    bool plink_format =false;
    int seed = 0;
    int thread = 1;
    double threshold = 0.0;
    int longIndex = 0;
    int opt = 0;
    opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while (opt != -1) {
        switch (opt)
        {
        case 'r': relate_name = optarg; break;
        case 'p': pheno_name = optarg; break;
        case 't':
            threshold = atof(optarg);
            if (threshold <= 0.0)
                fprintf(stderr,
                        "WARNING: Threshold = %f, will not filter samples\n",
                        threshold);
            break;
        case 's':
            try
            {
                int temp = misc::convert<int>(optarg);
                provide_seed = true;
                seed = temp;
            }
            catch (const std::runtime_error&)
            {
                fprintf(stderr, "Cannot parse the seed into number, will not "
                                "use the provided seed\n");
            }
            break;
        case 'o': out_name = optarg; break;
        case 'k': keep_name = optarg; break;
        case 'n':
            try
            {
                int temp = misc::convert<int>(optarg);
                thread = temp;
            }
            catch (const std::runtime_error&)
            {
                fprintf(stderr, "Cannot parse the thread into number, will "
                                "only use 1 thread\n");
            }
            if (thread <= 0) {
                fprintf(stderr, "Number of thread must be larger than 0. Will "
                                "use only 1 thread\n");
            }
            break;
        case 'i':
            id_1_col = optarg;
            plink_format = true;
            break;
        case 'I':
            id_2_col = optarg;
            plink_format = true;
            break;
        case 'f':
            f_col = optarg;
            plink_format = true;
            break;
        case 'h':
        case '?':
            usage();
            return 0;
        default:
            throw "Undefined operator, please use --help for more "
                  "information!";
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    if(use_plink_default){
        id_1_col = "IID1";
        id_2_col = "IID2";
        f_col="PI_HAT";
    }
    std::ostream* fp = &std::cout;
    std::ofstream fout;
    if (!out_name.empty()) {
        fout.open(out_name.c_str());
        if (!fout.is_open()) {
            fprintf(stderr, "ERROR: Cannot open output file to write: %s\n",
                    out_name.c_str());
            exit(-1);
        }
        fp = &fout;
    }
    std::string line;
    std::unordered_map<std::string, double> phenotype;
    std::unordered_set<std::string> include_samples;
    bool keep_samples = false;

    if(!keep_name.empty()){
        keep_samples = true;
        std::ifstream sample_file;
        sample_file.open(keep_name.c_str());
        if(!sample_file.is_open()){
            fprintf(stderr, "ERROR: Cannot open sample ID file from --keep:"
                            " %s\n", keep_name.c_str());
            exit(-1);
        }
        std::vector<std::string> token;
        while(getline(sample_file, line)){
            misc::trim(line);
            if(line.empty()) continue;
            token = misc::split(line);
            // we only use the first column of the file, not using both FID and IID
            include_samples.insert(token[0]);
        }
        sample_file.close();
    }
    if (!pheno_name.empty()) {
        std::ifstream pheno_file;
        pheno_file.open(pheno_name.c_str());
        if (!pheno_file.is_open()) {
            fprintf(stderr, "ERROR: Cannot open phenotype file %s\n",
                    pheno_name.c_str());
            exit(-1);
        }
        std::vector<std::string> token;
        while (getline(pheno_file, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            if (token.size() < 2) {
                fprintf(stderr, "ERROR: Phenotype file format incorrect! "
                                "Require at least 2 columns\n");
                exit(-1);
            }
            try
            {
                double factor = 0.0;
                if (token[1].compare("NA") == 0 || token[1].compare("na") == 0)
                    factor = -9;
                else
                    factor = misc::convert<double>(token[1]);
                if (phenotype.find(token[0]) != phenotype.end())
                    fprintf(stderr, "WARNING: Duplicated sample id: %s\n",
                            token[0].c_str());
                phenotype[token[0]] = factor;
            }
            catch (const std::runtime_error&)
            {
                // if it is not a number, we assume it is NA
                phenotype[token[0]] = -9;
            }
        }
        pheno_file.close();
    }

    std::random_device::result_type cur_seed = std::random_device()();
    if (provide_seed)
        cur_seed = static_cast<std::random_device::result_type>(seed);
    fprintf(stderr, "Seed used: %d\n", cur_seed);
    std::mt19937 rand_gen(cur_seed);
    // Read relationship file

    std::vector<Sample*> sample_list;
    if(!plink_format){
        sample_list=kin3col(relate_name, include_samples, phenotype,
                keep_samples, threshold, rand_gen);
    }else{
        sample_list=plink_format_process(relate_name, include_samples,
                                         phenotype, keep_samples, threshold, rand_gen,
                                         id_1_col, id_2_col, f_col);
    }
    // Update phenotype informations
    std::sort(sample_list.begin(), sample_list.end(), Sample::compare_sample);
    for (auto&& sample : sample_list) {
        if (sample->removed()) continue;
        sample->remove(*fp);
    }
    std::for_each(sample_list.begin(), sample_list.end(),
                  delete_pointed_to<Sample>);
    return 0;
}
