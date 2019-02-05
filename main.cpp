//
//  main.cpp
//  greedRelate
//  Greedy approach of keeping samples
//  Created by Shing Wan Choi on 18/08/2016.
//  Copyright © 2016 Shing Wan Choi. All rights reserved.
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
            if (!relative->removed()) {
                if (relative->m_occur > m_occur) {
                    relative->remove(os);
                }
                else if (relative->m_occur == m_occur)
                {
                    if (relative->m_phenotype < m_phenotype) {
                        relative->remove(os);
                    }
                    else if (misc::logically_equal(relative->m_phenotype,
                                                   m_phenotype))
                    {
                        if (relative->m_rand_number < m_rand_number) {
                            relative->remove(os);
                        }
                        else if (misc::logically_equal(relative->m_rand_number,
                                                       m_rand_number))
                        {
                            if (relative->m_name < m_name) {
                                relative->remove(os);
                            }
                        }
                    }
                }
            }
            if (m_occur <= 0 || m_removed) return 0;
        }
        // if we reach here, it means this sample still need to be removed.
        os << m_name << "\t" << m_occur << std::endl;
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
    fprintf(stderr, " v1.1.3 ( 2018-02-04 )\n");
    fprintf(stderr, " ==============================\n");
    fprintf(stderr, " This programme will try to minize the number of samples "
                    "that need to\n");
    fprintf(stderr, " be removed due to relatedness\n");
    fprintf(stderr,
            " Samples that should be removed are output to the STDOUT\n ");
    fprintf(stderr, " \n");
    fprintf(stderr, " Usage: GreedyRelate [options] -r <relatedness>\n");
    fprintf(stderr, "       -r | --relate      Relationship file (Required)\n");
    fprintf(stderr, "       -p | --pheno       Phenotype file\n");
    fprintf(stderr, "       -t | --threshold   Relatedness Threshold\n");
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


int main(int argc, char* argv[])
{
    if (argc <= 1) {
        usage();
        exit(0);
    }
    static const char* optString = "r:p:t:s:n:o:h?";
    static const struct option longOpts[] = {
        {"relate", required_argument, nullptr, 'r'},
        {"pheno", required_argument, nullptr, 'p'},
        {"threshold", required_argument, nullptr, 't'},
        {"thread", required_argument, nullptr, 'n'},
        {"seed", required_argument, nullptr, 't'},
        {"out", required_argument, nullptr, 'o'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};
    std::string relate_name = "";
    std::string pheno_name = "";
    std::string out_name = "";

    bool provide_seed = false;
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
        case 'h':
        case '?':
            usage();
            return 0;
            break;
        default:
            throw "Undefined operator, please use --help for more "
                  "information!";
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
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
    if (!pheno_name.empty()) {
        std::ifstream pheno_file;
        pheno_file.open(pheno_name.c_str());
        if (!pheno_file.is_open()) {
            fprintf(stderr, "ERROR: Cannot open phenotype file %s\n",
                    pheno_name.c_str());
            exit(-1);
        }
        // Assumine phenotype has header?
        getline(pheno_file, line);
        while (getline(pheno_file, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            std::vector<std::string> token = misc::split(line);
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
                fprintf(stderr, "ERROR: Undefined factor number\n");
            }
        }
        pheno_file.close();
    }

    std::random_device::result_type cur_seed = std::random_device()();
    ;
    if (provide_seed)
        cur_seed = static_cast<std::random_device::result_type>(seed);
    fprintf(stderr, "Seed used: %d\n", cur_seed);
    std::mt19937 rand_gen(cur_seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    // Read relationship file
    std::ifstream relate;
    relate.open(relate_name.c_str());
    if (!relate.is_open()) {
        fprintf(stderr, "ERROR: Cannot open relationship file %s\n",
                relate_name.c_str());
        exit(-1);
    }
    std::vector<Sample*> sample_list;
    std::unordered_map<std::string, size_t> sample_index;
    std::unordered_map<size_t, size_t> direction; // First size_t = pair, second
                                                  // size_t = index of the
                                                  // neighbour

    // Assume there is a header
    // Or we can add in
    fprintf(stderr,
            "Assuming there is a header file for the relatedness file\n");
    getline(relate, line);
    std::vector<std::string> token;
    size_t sample_idx = 0;
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
        double pheno = -9;
        if (phenotype.find(id) != phenotype.end()) pheno = phenotype[id];
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
    // Update phenotype informations
    sample_index.clear();
    std::sort(sample_list.begin(), sample_list.end(), Sample::compare_sample);
    for (auto&& sample : sample_list) {
        if (sample->removed()) continue;
        sample->remove(*fp);
    }
    std::for_each(sample_list.begin(), sample_list.end(),
                  delete_pointed_to<Sample>);
    return 0;
}
