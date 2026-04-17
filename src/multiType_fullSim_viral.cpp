//
//  multiType_FullSim.cpp
//
//
//  Created by Brandon Simony on 3/13/24.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <chrono>
#include <random>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <cassert>
#include <iomanip>
#include "gsl/gsl_randist.h"


//-------------------------//
// RV Generators using GSL //
//-------------------------//

//For initiating and seeding gsl random number generator.
class Gsl_rng_wrapper
{
    gsl_rng* r;
    public:
        Gsl_rng_wrapper()
        {
            std::random_device rng_dev; //To generate a safe seed.
            long seed = time(NULL)*rng_dev();
            const gsl_rng_type* rng_type = gsl_rng_default;
            r = gsl_rng_alloc(rng_type);
            gsl_rng_set(r, seed);
        }
        ~Gsl_rng_wrapper() { gsl_rng_free(r); }
        gsl_rng* get_r() { return r; }
};

// Long chunked Binomial RV using gsl.
long long draw_binom_gsl(long long N, double prob)
{
    static Gsl_rng_wrapper rng_w;
        static gsl_rng* r = rng_w.get_r();

        const unsigned int chunk_max = std::numeric_limits<unsigned int>::max();
        long long remaining = N;
        long long total_successes = 0;

        while (remaining > 0) {
            unsigned int chunk = (remaining > chunk_max) ? chunk_max : static_cast<unsigned int>(remaining);
            total_successes += gsl_ran_binomial(r, prob, chunk);
            remaining -= chunk;
        }

        return total_successes;
}


//Negative Binomial RV using gsl. -- parameterized using mean-dispersion as a generalization
// k = 1 is equivalent to a geometric rv and approaches a poisson rv as k -> infinity
// this process can equivalently be expressed as a poisson rv with a gamma distributed rate
// with shape parameter k and scale parameter R/k using the below parameterization
int draw_neg_binom_gsl(double R, double k)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    
    double p = 1 / (1 + (R / k));
    return gsl_ran_negative_binomial(r, p, k);
}

//Long chunked Multinomial RV using gsl.
std::vector<long long> draw_multinom_gsl(long long N, const std::vector<double>& prob)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    
    size_t K = prob.size();
    std::vector<long long> result(K, 0);
    const unsigned int chunk_max = std::numeric_limits<unsigned int>::max();

    long long remaining = N;
    while (remaining > 0) {
            unsigned int chunk = (remaining > chunk_max)
                               ? chunk_max
                               : static_cast<unsigned int>(remaining);

            std::vector<unsigned int> counts(K, 0);
            gsl_ran_multinomial(r, K, chunk, prob.data(), counts.data());

            // accumulate into long long result
            for (size_t i = 0; i < K; i++) {
                result[i] += counts[i];
            }

            remaining -= chunk;
        }

    return result;
}


//----------------------//
// Core model structure //
//----------------------//

// struct to hold parameter value inputs
struct Parameters
{
    int n_drugs = 1;
    int n_classes = 2;
    
    std::vector<long long> N_0; // initial population vector
    
    std::vector<double> b_vec; // type-specific birth rate
    std::vector<double> d_vec; // type-specific death rate
    std::vector<double> mu_vec; // per-site beneficial mutation probabilities
    std::vector<double> k_vec; // drug efficacy factor
    double mu_del; // probability of a deleterious mutation during reproduction
    double fit_cost; // per mutation cumulative fitness cost
    int n_retry; // number of allowed retries when extinction occurs before treatment
    double R; // mean viral burst size
    double k; // dispersion parameter of burst size distribution
    
    
    int T_max = 10; // number of generations
    int T_maxTreat = 10; // number of generations
    long long N_max = 1e9; // critical population size
    long long N_maxTreat = 2e9; // critical population size with treatment
    long long fullRes_nonExtinction = 1e8; // separate threshold for fully resistant mutants --> almost assured treatment failiure
    int verbose = 0; // debugging statements
    int data_out = 0; // specify model type. 0: console; 1: to file
    int reps_keep = 0; // number of replicates to store in full, not just endpoints
    std::string batchname = "outSize_distn"; // output file name tag
    std::string par_file = "pars.txt"; // output file name tag
    
    std::vector<std::string> conf_v;

    //Constructor either reads parameters from standard input (if no argument is given),
    //or from file (argument 1(.
    Parameters(int argc, char* argv[])
    {
        conf_v.reserve(200);
        std::stringstream buffer;
        if(argc == 3) //Config file
        {
            std::ifstream f(argv[2]);
            if(f.is_open())
            {
                buffer << f.rdbuf();
            }
            else
            {
                std::cout << "Failed to read config file \"" << argv[2] << "\"" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            for(int i=2; i<argc; ++i)
            {
                buffer << argv[i] << std::endl;
            }
        }
        while(!buffer.eof())
        { // until end of the stream
            std::string line = "";
            std::stringstream line_ss;
            // First get a whole line from the config file
            std::getline(buffer, line);
            // Put that in a stringstream (required for getline) and use getline again
            // to only read up until the comment character (delimiter).
            line_ss << line;
            std::getline(line_ss, line, '#');
            // If there is whitespace between the value of interest and the comment character
            // this will be read as part of the value. Therefore strip it of whitespace first.
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line.size() != 0)
            {
                conf_v.push_back(line);
            }
        }

        // must be read in before evaluating length since the number of classes (and thus parameters) is variable based on the number of variable sites / drugs
        
        size_t conf_length = 16; //number of model input parameters
        if (conf_v.size()!=conf_length)
        {
            std::cout << "Expected configuration file with " << conf_length << " options, loaded file with "
                      << conf_v.size() << " lines." << std::endl;
            exit(EXIT_FAILURE);
        }

        n_drugs = std::stoi(conf_v[0]); // number of drugs::resistance sites
        n_classes = pow(2, n_drugs);
        
        T_max = std::stoi(conf_v[1]); // simulation length
        N_max = std::stoll(conf_v[2]); // critical population size
        T_maxTreat = std::stoi(conf_v[3]); // select step between output generations
        N_maxTreat = std::stoll(conf_v[4]); // critical population size
        fullRes_nonExtinction = std::stoll(conf_v[5]); // separate threshold for fully resistant mutants --> almost assured treatment failiure
        mu_del = std::stod(conf_v[13]); // probability of a deleterious mutation during reproduction
        fit_cost = std::stod(conf_v[11]); // per mutation fitness cost
        n_retry = std::stoi(conf_v[12]); // number of allowed retries when extinction occurs before treatment
        R = std::stod(conf_v[14]); // mean viral burst size
        k = std::stod(conf_v[15]); // dispersion parameter of burst size distribution
        
        // simulation output
        verbose = std::stoi(conf_v[6]); // verbose statements (always written to batchname_err.txt file)
        data_out = std::stoi(conf_v[7]); // output indicator
        reps_keep = std::stoi(conf_v[8]); // number of replicates to store in full, not just endpoints
        batchname = conf_v[9]; // output file name tag
        par_file = conf_v[10]; // "path/f_name.txt" for vectorized parameters
        
        //-------------------------------//
        // read in vectorized parameters //
        //-------------------------------//
        
        std::string file_path = "./simulation_files/" + par_file + ".txt";
        std::ifstream parms(file_path);
        if (parms.is_open()) {
            std::string line;
            int lineIndex = 0;

            while (std::getline(parms, line)) {
                if (line.empty()) continue; // skip blank lines

                std::stringstream iss(line);
                std::string tmp;

                switch (lineIndex) {
                    case 0: // N_j
                        while (std::getline(iss, tmp, ',')) {
                            N_0.push_back(std::stoll(tmp));
                        }
                        break;
                    case 1: // b
                        while (std::getline(iss, tmp, ',')) {
                            b_vec.push_back(std::stod(tmp));
                        }
                        break;
                    case 2: // d
                        while (std::getline(iss, tmp, ',')) {
                            d_vec.push_back(std::stod(tmp));
                        }
                        break;
                    case 3: // mu
                        while (std::getline(iss, tmp, ',')) {
                            mu_vec.push_back(std::stod(tmp));
                        }
                        break;
                    case 4: // k
                        while (std::getline(iss, tmp, ',')) {
                            k_vec.push_back(std::stod(tmp));
                        }
                        break;
                    default:
                        std::cerr << "Warning: extra line in parameter file ignored: "
                                  << line << std::endl;
                        break;
                }

                lineIndex++;
            }

        } else {
            std::cout << "Failed to open file: " << file_path << std::endl;
            exit(EXIT_FAILURE);
        }
        
        parms.close();
        
        //-------------------------//
        // validate vector lengths //
        //-------------------------//
        
        // helper function for checks and error message
        auto check_size = [](const auto& vec, size_t expected,
                             const std::string& name, const std::string& ref,
                             std::ostringstream& err) {
            if (vec.size() != expected) {
                err << "Size mismatch: " << name << " (" << vec.size()
                    << ") vs " << ref << " (" << expected << ")\n";
            }
        };

        std::ostringstream err;
        
        // check for size mismatch
        check_size(k_vec, n_drugs, "k_vec", "n_drugs", err);
        check_size(mu_vec, n_drugs, "mu_vec", "n_drugs", err);
        check_size(N_0, n_classes, "N_0", "n_classes", err);
        check_size(b_vec, n_classes, "b_vec", "n_classes", err);
        check_size(d_vec, n_classes, "d_vec", "n_classes", err);
        

        // print all mismatch errors and terminate
        if (!err.str().empty()) {
            throw std::runtime_error(err.str());
        }

    }
};


// structure to contain population size of all classes at each generation
struct treeData {
    
    std::vector<std::vector<long long>> popData_full; // counts[class]
    std::vector<long long> N_vec;
    std::vector<int> treat_indicator;
    
    
    int genCount; // total number of generations in entire simulation
    int extinct; // indicator if simulation went extinct prior to treatment -- throw out replicate if == 1
    int crit; // indicator if emergence of total resistance occurred at any point -- count fraction of simulations == 1
    int treatSuccess; // indicator if treatment resolved infection -- crit + treatSuccess == 0 implies an unresolved simulation
    int preMutant;
    int preTreatCrit; // indicator if mutants were probabilistically critical (> p.fullRes_nonExtinction) before starting treatment.
    
    // ----- Default constructor -----
    treeData() :
        genCount(0),
        extinct(0),
        crit(0),
        treatSuccess(0),
        preMutant(0),
        preTreatCrit(0)
    {}
    
    treeData(int n_classes, int n_gen) {
        
        popData_full.assign(n_gen, std::vector<long long>(n_classes, 0LL)); // initialize counts to 0
        N_vec.assign(n_gen, 0LL);
        treat_indicator.assign(n_gen, 0);
        
        genCount = 0;
        extinct = 0;
        crit = 0;
        treatSuccess = 0;
        preMutant = 0;
        preTreatCrit = 0;
        
    }
    
    // Accessor (clearer than raw counts[g][c])
    long long& at(int gen, int cls) { return popData_full[gen][cls]; }
    
};


// simulation code for discrete branching process
// this corresponds to the discrete time generating function
treeData branchingSim_byGen( std::ofstream &errOut,
                             const Parameters& p,
                             const std::vector<std::vector<double>>& mu_prob,
                             const std::vector<double>& birth_probs_pre,
                             const std::vector<double>& birth_probs_post ){
    // initialize extinction reset variables
    int max_retries = p.n_retry;
    int retry_count = 0;
    treeData data;
    std::vector<long long> current_gen;
    std::vector<long long> next_gen;
    
    while(retry_count < max_retries){
        
        data = treeData(p.n_classes, p.T_max+p.T_maxTreat+1); // total generations + initial state for all class variables. last vector entry it total population size N
        
        // initialize counts, data structures at T = 0
        current_gen = p.N_0;
        next_gen.assign(current_gen.size(), 0);
        data.popData_full[0] = current_gen;
        data.N_vec[0] = std::accumulate(current_gen.begin(), current_gen.end(), 0LL);
        data.treat_indicator[0] = 0;
        data.genCount = 0;
        
        // iterate up to T_max or N_max -- pre-treatment
        for (int T = 1; T <= p.T_max; T++) {
            
            // Loop over types
            for (int type = 0; type < p.n_classes; ++type) {
                if (current_gen[type] == 0) continue;
                
                long long total_virions = 0;
                
                // each infected cell produces random number of free virions
                for(int parent = 0; parent < current_gen[type]; parent++){
                    total_virions += draw_neg_binom_gsl(p.R, p.k);
                }
                
                // number of free virions with lethal mutations
                total_virions -= (long long) draw_binom_gsl(total_virions, p.mu_del);
                
                // number of free virions that successfully infect a cell
                total_virions = draw_binom_gsl(total_virions, birth_probs_pre[type]);
                
               
                
                // number of offspring of each type - conditioned on not having a lethal mutation
                std::vector<long long> offspring_vector = draw_multinom_gsl(total_virions, std::vector<double>(mu_prob[type].begin(), mu_prob[type].end()) );
                
                // Add offspring to next generation
                for (int j = 0; j < offspring_vector.size(); ++j) {
                    next_gen[j] += offspring_vector[j];
                }
                
            }
            
            // stop early if negative value occurs, most likely due to large population overflow
            // abort before new (negative) values are added and replaced.
            if ( std::any_of(next_gen.begin(), next_gen.end(), [](long long x){ return x < 0LL; }) ) {
                errOut << "Suspected overflow detected due to negative class population. Aborting replicate simulation at generation T = " << T << std::endl;
                errOut << "Population size before error: " << data.N_vec[T-1] << std::endl;
                
                errOut << "Printing current population values: " << std::endl;
                for(int i = 0; i < current_gen.size(); i++){
                    errOut << current_gen[i] << "; ";
                }
                errOut << std::endl;
                
                errOut << "Printing next generation population values: " << std::endl;
                for(int i = 0; i < next_gen.size(); i++){
                    errOut << next_gen[i] << "; ";
                }
                errOut << std::endl;
                
                //--RESIZE NEW DATA STRUCTS--//
                data.popData_full.resize(data.genCount+1);
                data.N_vec.resize(data.genCount+1);
                data.treat_indicator.resize(data.genCount+1);
                return(data);
            }
            
            
            // store current counts
            data.N_vec[T] = std::accumulate(next_gen.begin(), next_gen.end(), 0LL); // total population size
            data.popData_full[T] = next_gen;
            data.treat_indicator[T] = 0;
            data.genCount++;
            
            
            // retry simulation if extinct (i.e., next gen. all zero). data.extinct now marks replicates that exceeded max_retries
            if ( data.N_vec[T] == 0 ) {
                
                retry_count ++; // increment retry count
                errOut << "Extinction occurred at pre-treatment generation " << T << std::endl;
                
                // terminate simulation if retry attempts exceeds threshold
                if (retry_count >= max_retries) {
                    errOut << "Maximum retries reached. Aborting simulation.\n";
                    data.extinct = 1;
                    
                    //--RESIZE NEW DATA STRUCTS--//
                    data.popData_full.resize(data.genCount+1);
                    data.N_vec.resize(data.genCount+1);
                    data.treat_indicator.resize(data.genCount+1);
                    
                    // reset current/future state storage containers
                    current_gen = next_gen; // next_gen becomes current_gen to start next iteration
                    next_gen.assign(next_gen.size(), 0); // next_gen is empty to start next iteration
                    
                    return(data);
                }
                
                continue;
                
            }
            
            // early stop if process exceeds population threshold -- pass endpoint into future sim with treatment
            if ( data.N_vec[T] >= p.N_max) {
                errOut << "Critical population of size " << data.N_vec[T] << " at generation " << T << std::endl;
                
                // reset current/future state storage containers
                current_gen = next_gen; // next_gen becomes current_gen to start next iteration
                next_gen.assign(next_gen.size(), 0); // next_gen is empty to start next iteration
                break;
            }
            
            // alternative stop if fully resistant population is almost surely critical (i.e., ((b+d) / b)^N_fullRes ~= 1) -- i.e., patient is kil
            if ( p.fullRes_nonExtinction != 0LL && next_gen[p.n_classes-1] >= p.fullRes_nonExtinction ) {
                data.crit = 1;
                data.preTreatCrit = 1;
                errOut << "Critical resistant population of size " << next_gen[p.n_classes-1] << " at generation " << T << std::endl;
                
                //--RESIZE NEW DATA STRUCTS--//
                data.popData_full.resize(data.genCount+1);
                data.N_vec.resize(data.genCount+1);
                data.treat_indicator.resize(data.genCount+1);
                
                // reset current/future state storage containers
                current_gen = next_gen; // next_gen becomes current_gen to start next iteration
                next_gen.assign(next_gen.size(), 0); // next_gen is empty to start next iteration
                
                // Check for mutants existing at critical point
                if (!current_gen.empty() && current_gen.back() > 0LL) {
                    data.preMutant = 2;
                }else if (std::any_of(current_gen.begin() + 1, current_gen.end(), [](long long x){ return x > 0LL; })) {
                    data.preMutant = 1;
                }
                
                return(data);
            }
            
            // reset current/future state storage containers
            current_gen = next_gen; // next_gen becomes current_gen to start next iteration
            next_gen.assign(next_gen.size(), 0); // next_gen is empty to start next iteration
            
        }
        
        // Check for mutants existing before treatment
        if (!current_gen.empty() && current_gen.back() > 0LL) {
            data.preMutant = 2;
            errOut << "Full mutant detected in population before treatment onset: " << current_gen.back() << std::endl;
        }else if (std::any_of(current_gen.begin() + 1, current_gen.end(), [](long long x){ return x > 0LL; })) {
            data.preMutant = 1;
        }
        
        break; // break from while loop if no other break conditions are triggered before treatment
    }
    
    int treatGen_count = 0; // count for treatment generations in case of pre-threshold critical or extinct state -- also serves to truncate excess output data
    int gen_preTreat = data.genCount; // effectively number of replication cycles before treatment begins -- also baseline indexing for post-treatment simulation data
    
    // TREATMENT SIMULATION
    // iterate to extinction or supercritical resistance
    for (int T = 1; T <= p.T_maxTreat; T++) {
        treatGen_count++;
        
        // Loop over types
        for (int type = 0; type < p.n_classes; ++type) {
            if (current_gen[type] == 0) continue;
            
            long long total_virions = 0;
            
            // each infected cell produces random number of free virions
            for(int parent = 0; parent < current_gen[type]; parent++){
                total_virions += draw_neg_binom_gsl(p.R, p.k);
            }
            
            // number of free virions with lethal mutations
            total_virions -= (long long) draw_binom_gsl(total_virions, p.mu_del);
            
            // number of free virions that successfully infect a cell
            total_virions = draw_binom_gsl(total_virions, birth_probs_post[type]);
            
           
            
            // number of offspring of each type - conditioned on not having a lethal mutation
            std::vector<long long> offspring_vector = draw_multinom_gsl(total_virions, std::vector<double>(mu_prob[type].begin(), mu_prob[type].end()) );
            
            // Add offspring to next generation
            for (int j = 0; j < offspring_vector.size(); ++j) {
                next_gen[j] += offspring_vector[j];
            }
            
        }

        // stop early if negative value occurs, most likely due to large population overflow
        // abort before new (negative) values are added and replaced.
        if ( std::any_of(next_gen.begin(), next_gen.end(), [](long long x){ return x < 0LL; }) ) {
            errOut << "Suspected overflow detected due to negative class population. Aborting replicate simulation at treatment generation T = " << T << std::endl;
            errOut << "Population size before error: " << data.N_vec[gen_preTreat + T] << std::endl;
            
            errOut << "Printing current population values: " << std::endl;
            for(int i = 0; i < current_gen.size(); i++){
                errOut << current_gen[i] << "; ";
            }
            errOut << std::endl;
            
            errOut << "Printing next generation population values: " << std::endl;
            for(int i = 0; i < next_gen.size(); i++){
                errOut << next_gen[i] << "; ";
            }
            errOut << std::endl;
            
            //--RESIZE NEW DATA STRUCTS--//
            data.popData_full.resize(data.genCount+1);
            data.N_vec.resize(data.genCount+1);
            data.treat_indicator.resize(data.genCount+1);
            return(data);
        }
        
        
        // store current counts
        data.N_vec[gen_preTreat + T] = std::accumulate(next_gen.begin(), next_gen.end(), 0LL); // total population size
        data.popData_full[gen_preTreat + T] = next_gen;
        data.treat_indicator[gen_preTreat + T] = 1;
        data.genCount++;
        
        
        // early stop if extinct (i.e., next gen. all zero), save extinction values
        if ( data.N_vec[gen_preTreat+T] == 0 ) {
            errOut << "Post treatment extinction occurred at treatment generation " << T << std::endl;
            data.treatSuccess = 1;
            // reset current/future state storage containers
            current_gen = next_gen; // next_gen becomes current_gen to start next iteration
            next_gen.assign(next_gen.size(), 0); // next_gen is empty to start next iteration
            break;
        }
        
        // early stop if process exceeds population threshold
        if ( data.N_vec[gen_preTreat + T] >= p.N_maxTreat) {
            data.crit = 1;
            errOut << "Critical post-treatment population of size " << data.N_vec[gen_preTreat + T] << " at treatment generation " << T << std::endl;
            // reset current/future state storage containers
            current_gen = next_gen; // next_gen becomes current_gen to start next iteration
            next_gen.assign(next_gen.size(), 0); // next_gen is empty to start next iteration
            break;
        }
        
        // alternative break if fully resistant population is almost surely critical (i.e., ((b+d) / b)^N_fullRes ~= 1) -- i.e., patient is kil
        if ( p.fullRes_nonExtinction != 0LL && next_gen[p.n_classes-1] >= p.fullRes_nonExtinction ) {
            data.crit = 1;
            errOut << "Critical resistant population of size " << next_gen[p.n_classes-1] << " at treatment generation " << T << std::endl;
            // reset current/future state storage containers
            current_gen = next_gen; // next_gen becomes current_gen to start next iteration
            next_gen.assign(next_gen.size(), 0); // next_gen is empty to start next iteration
            break;
        }
        
        // reset current/future state storage containers
        current_gen = next_gen; // next_gen becomes current_gen to start next iteration
        next_gen.assign(next_gen.size(), 0); // next_gen is empty to start next iteration
       
    }
    
    if(treatGen_count == p.T_maxTreat && data.N_vec[gen_preTreat + treatGen_count] > 0){
        errOut << "Simulation finished without resolution. Final size: " << data.N_vec[gen_preTreat + treatGen_count] << " at treatment generation " << treatGen_count << std::endl;
    }
    
    //--RESIZE NEW DATA STRUCTS--//
    errOut << "resizing simulation data from " << data.N_vec.size() << " to " << data.genCount+1 << std::endl;
    data.popData_full.resize(data.genCount+1);
    data.N_vec.resize(data.genCount+1);
    data.treat_indicator.resize(data.genCount+1);
    
    return(data);
    
}


//------------------------//
// Misc. helper functions //
//------------------------//

// Function to generate interpretable class names
std::vector<std::string> generateClassNames(int D) {
    std::vector<std::string> names;
    int total = 1 << D; // 2^K classes

    for (int mask = 0; mask < total; ++mask) {
        if (mask == 0) {
            names.push_back("WT"); // Wild type (no resistance)
        } else {
            std::string label = "M";
            for (int i = 0; i < D; ++i) {
                if (mask & (1 << i)) {
                    label += std::to_string(i + 1); // e.g., M1, M12
                }
            }
            names.push_back(label);
        }
    }
    return names;
}


// use bit-mapped type to determine number of mutations per individual
int numMutations(int x) {
    int count = 0;
    while (x) {
        count += x & 1;   // add 1 if lowest bit is set
        x >>= 1;          // shift right
    }
    return count;
}


// Function to compute offspring probability matrix
// P_{i,j} = P(type i has 1 type j offspring)
std::vector<std::vector<double>> offspringMatrix(const std::vector<double>& mu) {
    int K = mu.size();         // number of sites
    int numTypes = 1 << K; // 2^K possible genotypes -- bit representation
    
    std::vector<std::vector<double>> offspring_mat(numTypes, std::vector<double>(numTypes, -1.0));

    // For each parent genotype
    for (int parent = 0; parent < numTypes; ++parent) {
        // For each possible child genotype
        for (int child = 0; child < numTypes; ++child) {
            double prob = 1.0;
            // Loop over each site to compute probability of matching parent->child
            for (int site = 0; site < K; ++site) {
                bool parentHas = parent & (1 << site);
                bool childHas  = child  & (1 << site);

                if (parentHas == childHas) {
                    // no mutation at this site
                    prob *= (1.0 - mu[site]);
                } else {
                    // mutation occurred at this site
                    prob *= mu[site];
                }
            }
            offspring_mat[parent][child] = prob;
        }
    }
    return offspring_mat;
}


// Function for birth / death probabilities for each type given treatment
// calculates individual per-type birth probability when treatment is employed
// requires parameterization for alternative functional drug actions
std::vector<double> birthProbs(const Parameters& p, int treat_indicator){
    
    // declare b and d variables for
    // usage in lambda helper functions
    double b; double d;
    double pmc = p.fit_cost; // birth rate cost per mutation
    
    //------------------------------------//
    //     helper lambda functions for    //
    // different birth/death augmentation //
    //------------------------------------//
    
    // multiplies death by factor of k_j if susceptible
    auto f = [&b, &d](double k) -> std::vector<double> {
        double b_star = b;
        double d_star = k * d;
        return {b_star, d_star};
    };
    
    // reduction of birth -- WIP
    auto g = [&b, &d](double k) -> std::vector<double> {
        double b_star = b;
        double d_star = d;
        return {b_star, d_star};
    };
    
    // reduction of birth from fitness costs by # mutations
    // n_mut = numMutations(cls);
    auto fit_cost = [&b, &d, &pmc](int n_mut) -> std::vector<double> {
        double b_star = b*pow((1-pmc), n_mut);
        double d_star = d;
        return {b_star, d_star};
    };
    
    
    //---------------------//
    // start function body //
    //---------------------//
    
    std::vector<double> birthProbs;
    birthProbs.assign(p.n_classes, 0.0);
    
    for(int cls = 0; cls < p.n_classes; cls++) {
        
        // base birth/death
        b = p.b_vec[cls];
        d = p.d_vec[cls];
        
        // class-specific birth and death rates -- modified by fitness cost of mutation
        int cls_mut = numMutations(cls); // number of mutations on type cls individual
        auto bd_cost = fit_cost(cls_mut);
        b = bd_cost[0];
        d = bd_cost[1];
        
        // vector of birth probs under drug j treatment
        std::vector<double> typej_probs(p.n_drugs, 0.0);

        for(int drugType = 0; drugType < p.n_drugs; drugType++){
            
            // check for resistance to current drug type
            if( (cls >> drugType) & 1 ){
                // resistant, so no change
                typej_probs[drugType] = b / (b + d);
            }else{
                
                // only apply drug effects during treatment phase
                if(treat_indicator == 1){
                    auto new_bd = f(p.k_vec[drugType]); // susceptible, so apply augmentation
                    typej_probs[drugType] = new_bd[0] / (new_bd[0] + new_bd[1]);
                }else{
                    typej_probs[drugType] = b / (b + d);
                }
                
            }
            
        } // end loop over drug types
        
        birthProbs[cls] = *std::min_element(typej_probs.begin(), typej_probs.end()); // keep minimum birth probability -- implies only the strongest drug applies mortality
        
    } // end loop over classes
    
    return birthProbs;
    
}


//------------------------------------------//
// Function main performs parameter read-in //
// and manages replicate results and output //
//------------------------------------------//

int main(int argc, char* argv[]){
    
    // initialize parameters read from command line input
    int n_reps = std::stoi(argv[1]); // simulation replicates
    Parameters p(argc, argv);
    
    
    if(p.verbose > 2){ std::cout << std::filesystem::current_path() << std::endl; }
    
    // open error file
    std::string err_fname = "./simulation_files/" + p.batchname + "_err.txt";
    std::ofstream errOut(err_fname);
    if (!errOut.is_open()) {
        std::cerr << "Unable to open error file: " << err_fname << std::endl;
        return EXIT_FAILURE;
    }
    std::cerr.rdbuf(errOut.rdbuf());  // Redirect cerr to the error file
    
    // open output file
    std::string fname = "./simulation_files/" + p.batchname + "_results.txt";
    std::ofstream simResults(fname);
    if (!simResults.is_open()) {
        std::cerr << "Unable to open error file: " << fname << std::endl;
        return EXIT_FAILURE;
    }
    
    // open full trajectories output file
    std::string full_fname = "./simulation_files/" + p.batchname + "_full_results.txt";
    std::ofstream simFullResults(full_fname);
    if (!simFullResults.is_open()) {
        std::cerr << "Unable to open error file: " << full_fname << std::endl;
        return EXIT_FAILURE;
    }
    
    // write current runtime to .err file
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    errOut << std::put_time(std::localtime(&now_time), "%Y-%m-%d %H:%M:%S") << std::endl;
    // std::cout << "Wrote to " << err_fname << std::endl;
    
    // print non-vector parameters
    if(p.verbose > 0){
        
        errOut << "###################################\n"; // end of row
        errOut << "Printing non-vector parameters \n"; // end of row
        errOut << "###################################\n"; // end of row
        errOut << "n_drugs: " << p.n_drugs << "\n";
        errOut << "n_classes: " << p.n_classes << "\n";
        errOut << "T_max: " << p.T_max << "\n";
        errOut << "N_max: " << p.N_max << "\n";
        errOut << "T_maxTreat: " << p.T_maxTreat << "\n";
        errOut << "N_maxTreat: " << p.N_maxTreat << "\n";
        errOut << "fullRes_nonExtinction: " << p.fullRes_nonExtinction << "\n";
        errOut << "mu_del: " << p.mu_del << "\n";
        errOut << "verbose: " << p.verbose << "\n";
        errOut << "data_out: " << p.data_out << "\n";
        errOut << "batchname: " << p.batchname << "\n";
        errOut << "par_dat: " << p.par_file << "\n" << "\n" << std::endl;

    }
   
    // print vectorized parameters
    if(p.verbose > 0){
        
        errOut << "#############################\n"; // end of row
        errOut << "Printing vector parameters \n"; // end of row
        errOut << "#############################\n"; // end of row
        
        errOut << "N_0: ";
        for (size_t i = 0; i < p.N_0.size(); ++i) { errOut << p.N_0[i] << "; "; }
        errOut << "\n";
        
        errOut << "b_vec: ";
        for (size_t i = 0; i < p.b_vec.size(); ++i) { errOut << p.b_vec[i] << "; "; }
        errOut << "\n";
        
        errOut << "d_vec: ";
        for (size_t i = 0; i < p.d_vec.size(); ++i) { errOut << p.d_vec[i] << "; "; }
        errOut << "\n";
        
        errOut << "mu_vec: ";
        for (size_t i = 0; i < p.mu_vec.size(); ++i) { errOut << p.mu_vec[i] << "; "; }
        errOut << "\n";
        
        errOut << "k_vec: ";
        for (size_t i = 0; i < p.k_vec.size(); ++i) { errOut << p.k_vec[i] << "; "; }
        errOut << "\n" << "\n" << std::endl;
        
    }
    
    
    //-----------------------------//
    // pre-simulation computations //
    //-----------------------------//
    
    // calculate by-type offspring probability distribution  matrix
    std::vector<std::vector<double>> prob_matrix = offspringMatrix(p.mu_vec);
    
        
    // troubleshoot/debug visual for offspring matrix
    if(p.verbose > 0){
        
        errOut << "#####################################\n"; // end of row
        errOut << "Printing offspring prob. matrix\n"; // end of row
        errOut << "#####################################\n"; // end of row
        for (size_t i = 0; i < prob_matrix.size(); ++i) {
            for (size_t j = 0; j < prob_matrix[i].size(); ++j) {
                errOut << std::setw(10)               // fix column width
                       << std::scientific << std::setprecision(3)        // 4 digits after decimal
                       << prob_matrix[i][j];
                if (j < prob_matrix[i].size() - 1) errOut << " ";
            }
            errOut << "\n";
        }
        errOut << "\n" << std::endl;
    }
    
    // calculate by type birth probability before and after treatment initiated
    std::vector<double> preTreat_birth_vec  = birthProbs( p, 0 );
    std::vector<double> postTreat_birth_vec = birthProbs( p, 1 );
        
    // troubleshoot/debug visual for preTreatment birth vector
    if(p.verbose > 0){
        
        errOut << "###########################################\n"; // end of row
        errOut << "Printing pre-treatment birth prob. vector  \n"; // end of row
        errOut << "###########################################\n"; // end of row
        for (size_t i = 0; i < preTreat_birth_vec.size(); ++i) {
            errOut << preTreat_birth_vec[i] << "; ";
        }
        errOut << "\n" << "\n" << std::endl; // end of row
    }
    // troubleshoot/debug visual for postTreatment birth vector
    if(p.verbose > 0){
        
        errOut << "###########################################\n"; // end of row
        errOut << "Printing post-treatment birth prob. vector \n"; // end of row
        errOut << "###########################################\n"; // end of row
        for (size_t i = 0; i < postTreat_birth_vec.size(); ++i) {
            errOut << postTreat_birth_vec[i] << "; ";
        }
        errOut << "\n" << "\n" << std::endl; // end of row
    }
    
    //------------------------------------//
    // generate dynamic labels for header //
    //------------------------------------//
    
    std::vector<std::string> class_names = generateClassNames(p.n_drugs); // generate all class names for n_drugs drugs
    if(p.verbose > 0){
        
        errOut << "#######################\n"; // end of row
        errOut << "Printing class names\n"; // end of row
        errOut << "#######################\n"; // end of row
        for (size_t i = 0; i < class_names.size(); ++i) {
            errOut << class_names[i] << "; ";
        }
        errOut << "\n" << "\n" << std::endl; // end of row
    }
    

    //----------------------------//
    // select model output format //
    //----------------------------//
    
    // error file start simulations header
    if(p.verbose > 0){
        
        errOut << "#######################\n"; // end of row
        errOut << "Running Simulations\n"; // end of row
        errOut << "#######################\n"; // end of row
    }
    
            
    //------------------------------------//
    // run model and print output to file //
    //------------------------------------//
    
    // column headers
    simResults << "rep;" << "N;" << "ext;" << "crit;" << "preTreatCrit;" << "treatSuccess;" << "preResistance";
    for (int cls = 0; cls < p.n_classes; cls++) {
        simResults << ";" << class_names[cls];
    }
    simResults << std::endl;
    
    
    // column headers
    simFullResults << "rep;" << "Gen;" << "treatStatus";
    for (int cls = 0; cls < p.n_classes; cls++) {
        simFullResults << ";" << class_names[cls];
    }
    simFullResults << ";N" << std::endl;
    
    
    // model results by replicate
    int reps_stored = 0;
    for(int i = 0; i < n_reps; i++){
        // run model
        if(p.verbose > 0){ errOut << "Rep. " << i+1 << "\n"; }
        treeData res_out = branchingSim_byGen( errOut, p, prob_matrix, preTreat_birth_vec, postTreat_birth_vec );
        
        // errOut << "vector size: " << res_out.treat_indicator.size() << " v.s. genCount value: " << res_out.genCount << std::endl;
        
        // write out select replicate trajectories in full if treatment was reached
        if(p.reps_keep > reps_stored && res_out.extinct != 1 ){
            reps_stored++; // increment stored rep count
            
            // write out full replicate trajectory
            for(int gen = 0; gen <= res_out.genCount; gen++){
                simFullResults << i+1 << ";" << gen << ";" << res_out.treat_indicator[gen];
                for (int cls = 0; cls < p.n_classes; cls++) {
                    simFullResults << ";" << res_out.at(gen, cls);
                }
                simFullResults << ";" << res_out.N_vec[gen] << "\n";
            }
        }
        
        
        // post treatment endpoint values
        simResults << i+1 << ";" << res_out.N_vec[res_out.genCount] << ";" << res_out.extinct << ";" << res_out.crit << ";" << res_out.preTreatCrit << ";" << res_out.treatSuccess << ";" << res_out.preMutant;
        for (int cls = 0; cls < p.n_classes; cls++) {
            simResults << ";" << res_out.at(res_out.genCount, cls);
        }
        simResults << "\n";
    }
    
    simResults.close(); // result file
    simFullResults.close(); // select full trajectory file
    errOut.close(); // error file

    return 0;
    
}
