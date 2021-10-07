/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   manutility.h
 * Author: ali_nayeem
 *
 * Created on August 11, 2019, 12:38 PM
 */

#ifndef MANUTILITY_H
#define MANUTILITY_H

#include <cstdlib>
//#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>

#include "msa.h"

using namespace std;

std::string GetStdoutFromCommand(std::string cmd)
{

    std::string data;
    FILE *stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    //cmd.append(" 2>&1");

    stream = popen(cmd.c_str(), "r");
    if (stream)
    {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL)
                data.append(buffer);
        pclose(stream);
    }
    return data;
}

double softsign(double x)
{
    return x / (1.0 + abs(x));
}

unsigned countDigit(unsigned num)
{
    unsigned count = 0;
    while (num != 0)
    {
        /* Increment digit count */
        count++;
        /* Remove last digit of 'num' */
        num /= 10;
    }
    return count;
}

vector<unsigned> calculateDigitCount(vector<double> objVec)
{
    vector<unsigned> digVec;
    for (size_t i = 0; i < objVec.size(); i++)
    {
        digVec.push_back(countDigit((unsigned)objVec[i]));
    }
    return digVec;
}

unsigned calculateMaxDigitCount(vector<double> objVec)
{
    //vector<unsigned> digVec = calculateDigitCount(objVec);
    return countDigit((unsigned)*max_element(objVec.begin(), objVec.end()));
}
vector<double> calculateSimgSimngScore(const MSA &msa) //both maximization
{
    unsigned uColCount = msa.GetColCount();
    unsigned uRowCount = msa.GetSeqCount();
    unsigned gapColCount = 0, nonGapColCount = 0;
    double simg = 0, simng = 0;
    for (unsigned i = 0; i < uColCount; i++)
    {
        bool gapColumn = false;
        std::unordered_map<char, unsigned> columnMap;
        for (unsigned j = 0; j < uRowCount; j++)
        {
            if (msa.IsGap(j, i))
            {
                gapColumn = true;
            }
            else
            {
                columnMap[msa.GetChar(j, i)]++;
            }
        }
        auto x = std::max_element(columnMap.begin(), columnMap.end(), [](const std::pair<char, unsigned> &p1, const std::pair<char, unsigned> &p2) { return p1.second < p2.second; });
        unsigned max = x->second;
        if (gapColumn)
        {
            gapColCount++;
            simg += 1.0 * max / uRowCount;
        }
        else
        {
            nonGapColCount++;
            simng += 1.0 * max / uRowCount;
        }
    }
    vector<double> ret; //two values to be returned, STL pair could be used as well
    ret.push_back(simg);
    ret.push_back(simng);
    return ret;
    // simg = simg / gapColCount; //to scale below 1.0
    // simng = simng / nonGapColCount; //to scale below 1.0
    // SCORE result = simg * simgWeight + simng * simngWeight;
    //printf("c++=%lf, %lf, %lf\n", simg, simng, result);
    //return -1 * result; //possibly a bug from MAN: both simg and simng are maximizing, seems -1 not necessary
}

double calculateGapScore(const MSA &msa) //minimization
{
    unsigned uColCount = msa.GetColCount();
    unsigned uRowCount = msa.GetSeqCount();

    double gapScoreSum = 0.0;

    for (unsigned i = 0; i < uRowCount; i++)
    {
        unsigned gapCount = 0;
        unsigned nonGapCount = 0;
        for (unsigned j = 0; j < uColCount; j++)
        {
            if (msa.IsGap(i, j))
                gapCount++;
            else
                nonGapCount++;
        }
        gapScoreSum += gapCount;
        //gapScoreSum += 1.0 * gapCount / nonGapCount; //divide by nonGapCount is used to scale below 1.0
    }

    return gapScoreSum;
    //return -gapWeight * ( 1.0 * gapScoreSum / uRowCount ) ; //averaging used to scale below 1.0
}

double aggregationFunction(vector<double> objVector, vector<float> & weightVector, int functionType_ = 0)
{
    double score;
    if (functionType_ == 0)
    {
        score = inner_product(objVector.begin(), objVector.end(), weightVector.begin(), 0.0); 
    }
    else
    {
        //printf("In TCHEB");
        std::transform(objVector.begin(), objVector.end(), objVector.begin(), [&](double x){return(1.1 - x);});
        vector<double> result(objVector.size());
        std::transform(objVector.begin(), objVector.end(), weightVector.begin(), result.begin(), std::multiplies<double>());
        score = -1.0 * (*max_element(result.begin(), result.end())); //bcoz Tchebycheff gives a minimization score, but MUSCLE tries to maximize
        // double maxFun = -1.0e+30;

        // for (int n = 0; n < objVector.size(); n++)
        // {
        //     double objValue = objVector[n];

        //     double feval;
        //     if (weightVector[n] == 0.0)
        //     {
        //         feval = 0.00001 * objValue;
        //     }
        //     else
        //     {
        //         feval = objValue * weightVector[n];
        //     }
        //     if (feval > maxFun)
        //     {
        //         maxFun = feval;
        //     }
        // } // for

        //fitness = maxFun;
    } // if
    return score;
}
#endif /* MANUTILITY_H */
