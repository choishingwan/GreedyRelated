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


class Sample
{
public:
    // First occurance = itself, therefore occurance = 0
    Sample(std::string name, double pheno, double rand)
        : m_name(name), phenotype(pheno), rand_number(rand)
    {
        occur = 0;
    }
    // add a relative to the current sample
    int add(Sample* related)
    {
        // not sure why
        // TODO: See if we can skip the step of pushing nullptr
        relatives.push_back(nullptr);
        relatives.back() = related;
        m_removed = false;
        occur++;
        return occur;
    }
    int remove()
    {
        // we are going to start cleaning up these relatives.
        // good thing is that the size of this vector will small
        // check if any relatives have more than us. If not, then proceed,
        // otherwise, return
        if (m_name == "3003339" || m_name == "5128261" || m_name == "1304775"
            || m_name == "5000708" || m_name == "2975133")
        {
            std::cerr << "Removing individual: " << m_name << "\t" << occur
                      << std::endl;
        }
        for (auto&& relative : relatives) {
            if (!relative->removed()
                && (relative->occur > occur
                    || (relative->occur == occur
                        && relative->phenotype < phenotype)
                    || (relative->occur == occur
                        && relative->rand_number > rand_number)
                    || (relative->occur == occur
                        && relative->rand_number == rand_number
                        && relative->m_name > m_name)))
            {

                if (m_name == "3003339" || m_name == "5128261"
                    || m_name == "1304775" || m_name == "5000708"
                    || m_name == "2975133")
                {
                    std::cerr << "Relative Save: " << relative->m_name << "\t"
                              << relative->occur << "\t"
                              << relative->rand_number << "\t" << rand_number
                              << std::endl;
                }
                return 0;
            }
        }
        std::cout << m_name << "\t" << occur << std::endl;
        occur = -1;
        m_removed = true;
        for (auto&& relative : relatives) {
            relative->occur--;
            if (relative->occur <= 0) relative->m_removed = true;
            if (relative->m_name == "3003339" || relative->m_name == "5128261"
                || relative->m_name == "1304775"
                || relative->m_name == "5000708"
                || relative->m_name == "2975133")
            {
                std::cerr << "Relative reduce: " << m_name << "\t" << occur
                          << "\t" << relative->m_name << "\t" << relative->occur
                          << "\t" << relative->rand_number << "\t"
                          << rand_number << std::endl;
            }
        }
        std::sort(relatives.begin(), relatives.end(), Sample::compare_sample);
        for (auto&& relative : relatives) {
            if (!relative->removed() && relative->occur > 0) {
                if (relative->m_name == "3003339"
                    || relative->m_name == "5128261"
                    || relative->m_name == "1304775"
                    || relative->m_name == "5000708"
                    || relative->m_name == "2975133")
                {
                    std::cerr << "Start remove relative" << std::endl;
                }
                relative->remove();
            }
        }
        return occur;
    }

    std::string debug() const
    {
        std::string ocur;           // string which will contain the result
        std::ostringstream convert; // stream used for the conversion
        convert << occur; // insert the textual representation of 'Number' in
                          // the characters in the stream
        ocur = convert.str();
        return m_name + " " + ocur;
    }

    static bool compare_sample(Sample const* a, Sample const* b)
    {
        if (a->occur == b->occur) {
            if (a->phenotype == b->phenotype) {
                if (a->rand_number == b->rand_number) {
                    return a->m_name > b->m_name;
                }
                else
                    return a->rand_number > b->rand_number;
            }
            else
                return a->phenotype < b->phenotype;
        }
        else
            return a->occur > b->occur;
    }

    int get_occur() const { return occur; }
    std::string get_name() const { return m_name; }
    bool removed() const
    {
        // add additional condition where the current sample is removed
        if (occur == 0) return true;
        return m_removed;
    }

private:
    std::vector<Sample*> relatives;
    std::string m_name;
    double phenotype;
    double rand_number; // This is a random number to solve the tie
    int occur;
    bool m_removed = true; // so any samples with invalid pairs will be ignoreds
};

void usage()
{
    fprintf(stderr, " GreedyRelate\n");
    fprintf(stderr, " Sam Choi\n");
    fprintf(stderr, " v1.1.1 ( 2018-01-23 )\n");
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

int main(int argc, char* argv[])
{
    if (argc <= 1) {
        usage();
        exit(0);
    }
    static const char* optString = "r:p:t:s:n:h?";
    static const struct option longOpts[] = {
        {"relate", required_argument, nullptr, 'r'},
        {"pheno", required_argument, nullptr, 'p'},
        {"threshold", required_argument, nullptr, 't'},
        {"thread", required_argument, nullptr, 'n'},
        {"seed", required_argument, nullptr, 't'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, 0, 0}};
    std::string relate_name = "";
    std::string pheno_name = "";
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
            catch (const std::runtime_error& error)
            {
                fprintf(stderr, "ERROR: Undefined factor number\n");
            }
        }
        pheno_file.close();
    }

    std::random_device::result_type cur_seed = std::random_device()();
    ;
    if (provide_seed) cur_seed = seed;
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
        catch (const std::runtime_error& error)
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
        sample->remove();
    }
    return 0;
}
