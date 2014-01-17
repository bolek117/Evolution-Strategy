#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>
#include <sstream>

#define uint unsigned int
#define PI 3.14159265359f

using namespace std;

float randn();
vector < vector < float > > loadData(string filename);
float f(float a, float b, float c, float x);
bool sortByColumn(const vector < float > &v1, const vector < float > &v2);
bool saveData(bool isHeader, uint mi=0, uint lambda=0, uint k=0, string strategy="", float eps=0.0f, float a=0.0f, float b=0.0f, float c=0.0f, uint time=0, unsigned int iterationsToFind=0);

enum Types
{
    CHROMOSOME = 0,
    STDDEV
};

enum Strategies
{
    MI_COMMA_LAMBDA = 0,
    MI_PLUS_LAMBDA
};

int main()
{
    vector < uint > miV, lambdaV, kV, strategyV;
    miV.push_back(10);
    miV.push_back(25);
    miV.push_back(50);
    miV.push_back(100);
    miV.push_back(250);

    lambdaV.push_back(2);
    lambdaV.push_back(10);
    lambdaV.push_back(50);

    kV.push_back(100);
    kV.push_back(500);
    kV.push_back(1000);
    kV.push_back(2500);

    strategyV.push_back(MI_PLUS_LAMBDA);
    strategyV.push_back(MI_COMMA_LAMBDA);


    float eps = 0.215f,
        max = 10.0f;

    const uint n=3; // DO NOT CHANGE!

    uint mi=100, lambda=2, K=500, strategy=MI_PLUS_LAMBDA;
    string strategyStr;
    saveData(true, 0, 0, 0, "", 0);
    uint startStart = time(NULL);
    for (uint iStrategy=0;iStrategy<strategyV.size();++iStrategy)
    {
        strategy = strategyV[iStrategy];
        if (strategy == MI_COMMA_LAMBDA)
            strategyStr = "mi_comma_lambda";
        else
            strategyStr = "mi_plus_lambda";

        for (uint iK=0;iK<kV.size();++iK)
        {
            K = kV[iK];

            for (uint iMi=0;iMi<miV.size();++iMi)
            {
                mi = miV[iMi];
                for (uint iLambda=0;iLambda<lambdaV.size();++iLambda)
                {
                    lambda = lambdaV[iLambda];

                    cout << "Strategy: " << strategyStr << ", maxIterations: " << K << ", parents: " << mi << ", childs: " << lambda << "\n";

                    // init some variables for performance purposes (variables should be not realocated)
                    srand((uint)time(NULL));
                    vector < float > bestSolution;
                    vector < vector < float > > data;

                    vector < float > genes;
                    genes.assign(n , 0.0f);	// a, b, c

                    vector < vector < float > > chromosomes;
                    chromosomes.assign(2, genes);	// chromosome & standard deviation

                    vector < float> pair;
                    pair.assign(2, 0.0f);

                    vector < float > lastChange;
                    lastChange = pair;

                    vector < vector < float > > means;	// mean of population for quality evaluation
                    vector < vector < float > > sorted;	// for sorted means

                    vector < vector < vector < float > > > population;
                    population.assign(mi, chromosomes);	// all chromosome & stdDev pairs

                    vector < vector < vector < float > > > nextGeneration;
                    nextGeneration.assign(mi*lambda, chromosomes);	// vector for childrens

                    vector < vector < vector < float > > > finalPopulation;
                    finalPopulation = population;

                    // load data from input file and handle file reading errors
                    data = loadData("data.dat");
                    if (data.size() == 0)
                    {
                        cout << "[ERROR] File \"data.dat\" not found or not readable. Terminating...";
                        return -1;
                    }

                    time_t start = time(NULL);

                    // initial population generation (parent -> elementType(chromosome/stdDev) -> chromosome -> gene)
                    for (uint parent=0;parent<mi;++parent)
                    {
                        for (uint elementType=0;elementType<2;++elementType)
                        {
                            for (uint gene=0;gene<n;++gene)
                            {
                                float rand = randn();

                                if (elementType == STDDEV)
                                    rand = abs(rand);

                                if (abs(rand) <= max)
                                    genes.at(gene) = rand;
                                else
                                    --gene;
                            }

                            chromosomes.at(elementType) = genes;
                        }

                        population.at(parent) = chromosomes;
                    }

                    // evaluate TAU's
                    float tau1, tau2;
                    tau1 = 1.0f/sqrt(2.0f * (float)n);
                    tau2 = 1.0f/sqrt(2.0f * sqrt((float)n));

                    // evaluate next generation
                    float bestQuality = 9.99e30f;	// initialize best quality variable with some huge variable (near max of float)

                    // initialize help variables
                    float nextDev, nextGene, actualGene, actualDev, sum, x, y, fVal, mean, expRandnTau1;
                    uint index;

                    for (uint k=0;k<K && bestQuality>eps;++k)
                    {
                        // for every parent
                        index = 0;
                        for (uint parent=0;parent<mi;++parent)
                        {
                            // make 'lambda' of childs
                            for (uint child=0;child<lambda;++child)
                            {
                                expRandnTau1 = exp(randn()*tau1);

                                // for every gene
                                for (uint gene=0;gene<n;++gene)
                                {
                                    actualGene = population.at(parent).at(CHROMOSOME).at(gene);
                                    actualDev = population.at(parent).at(STDDEV).at(gene);

                                    nextDev = actualDev * expRandnTau1 * exp(randn()*tau2);
                                    nextGene = actualGene + randn() * nextDev;

                                    chromosomes.at(CHROMOSOME).at(gene) = nextGene;
                                    chromosomes.at(STDDEV).at(gene) = nextDev;
                                }

                                nextGeneration.at(index) = chromosomes;
                                ++index;
                            }
                        }

                        if (strategy == MI_PLUS_LAMBDA)
                            population.insert(population.end(), nextGeneration.begin(), nextGeneration.end());
                        else
                            population = nextGeneration;

                        // for every specimen evaluate bi-quadratic mean
                        genes.assign(n, 0.0f);
                        means.assign(population.size(), pair);

                        for (uint specimen=0;specimen<population.size();++specimen)
                        {
                            // calculate sum for bi-quadratic mean
                            genes = population.at(specimen).at(CHROMOSOME);
                            sum = 0;

                            for (uint i=0;i<data.size();++i)
                            {
                                x = data.at(i).at(0);
                                y = data.at(i).at(1);
                                fVal = y - f(genes.at(0), genes.at(1), genes.at(2), x);
                                sum += powf(fVal, 2.0f);
                            }

                            mean = sum/data.size();
                            means.at(specimen).at(0) = mean;
                            means.at(specimen).at(1) = (float)specimen;
                        }

                        // sort `means` for selecting the best
                        sort(means.begin(), means.end(), sortByColumn);

                        bestQuality = means.at(0).at(0);
                        bestSolution = population.at((uint)means.at(0).at(1)).at(CHROMOSOME);

                        if (bestQuality != lastChange.at(0))
                        {
                            lastChange.at(0) = bestQuality;
                            lastChange.at(1) = (float)k;
                        }

                        for (uint i=0;i<mi;++i)
                            finalPopulation.at(i) = population.at((uint)means.at(i).at(1));

                        population = finalPopulation;

                        if (k%(K/50) == 0 || k == K-1)
                        {
                            cout << "|";
                            cout.flush();
                            //cout << "i=" << k << ", q=" << bestQuality << ", a=" << bestSolution.at(0) << ", b=" << bestSolution.at(1) << ", c=" << bestSolution.at(2) <<"\n";
                        }
                    }

                    cout << "\n";
                    time_t end = time(NULL);

                    cout << "Total time: " << end-start << "s\n"
                        << "Best solution in itearation: " << (uint)lastChange.at(1) << "\n";

                    saveData(false, mi, lambda, K, strategyStr, eps, bestSolution.at(0), bestSolution.at(1), bestSolution.at(2), end-start, (uint)lastChange.at(1));

                    cout << "\n-----------------------------------------\n\n";
                }
            }
        }
    }

    uint stopStop = time(NULL);
    cout << "Total time: " << stopStop - startStart << "s\n";
    return 0;
}

float randn()
{
    float v1,v2,s;

    do
    {
        v1 = 2.0f * ((float) rand()/RAND_MAX) - 1;
        v2 = 2.0f * ((float) rand()/RAND_MAX) - 1;

        s = powf(v1, 2.0f) + powf(v2, 2.0f);
    }
    while ( s >= 1.0f );

    if (s == 0.0f)
        return 0.0f;
    else
        return (v1*sqrt(-2.0f * log(s) / s));
}

float f(float a, float b, float c, float x)
{
    return a * (pow(x, 2.0f) - b * cos(c*PI*x));
}

vector < vector < float > > loadData(string filename)
{
    ifstream data;
    data.open(filename.c_str());
    char line[64];
    string tmpS;
    vector < vector < float > > d;
    if (data.is_open())
    {
        while (data.getline(line, 64))
        {
            vector < float > tmp0;

            tmpS = string(line);

            // delete spaces from begining
            for (uint i=0;i<tmpS.size();++i)
            {
                if (tmpS[0] == 0x20)
                    tmpS = tmpS.substr(1);
                else
                    break;
            }

            tmp0.push_back((float)atof(tmpS.c_str()));
            tmp0.push_back((float)atof(tmpS.substr(tmpS.find_first_of(' ')).c_str()));

            d.push_back(tmp0);
        }

        return d;
    }
    else
    {
        vector < vector < float > > emptyRes;
        return emptyRes;
    }
}

bool sortByColumn(const vector < float > &v1, const vector < float > &v2)
{
    return abs(v1[0]) < abs(v2[0]);
}

// TODO save data
bool saveData(bool isHeader, uint mi, uint lambda, uint k, string strategy, float eps, float a, float b, float c, uint time, uint iterationsToFind)
{
    ofstream file;
    stringstream out;

    if (!isHeader)
    {
        file.open("output.csv", ofstream::app);
        out << mi << ";"
            << lambda << ";"
            << k << ";"
            << strategy << ";"
            << eps << ";"
            << a << ";"
            << b << ";"
            << c << ";"
            << time << ";"
            << iterationsToFind << "\n";
    }
    else
    {
        file.open("output.csv", ofstream::trunc);
        out << "parents;childs;k;strategy;eps;a;b;c;time;iterationsToFind\n";
    }

    if (file.good())
    {
        file << out.str();
        return true;
    }
    else
        return false;
}
