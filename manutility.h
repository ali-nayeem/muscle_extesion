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
#include<unordered_map>
#include <algorithm>
#include <vector>

#include "msa.h"

using namespace std;

std::string GetStdoutFromCommand(std::string cmd) {

    std::string data;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    //cmd.append(" 2>&1");

    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
        pclose(stream);
    }
    return data;
}

double softsign(double x)
{
    return x / (1.0 + abs(x));
}

vector<double> calculateSimgSimngScore(const MSA &msa) //both maximization
{
    unsigned uColCount = msa.GetColCount();
    unsigned uRowCount = msa.GetSeqCount();
    unsigned gapColCount = 0, nonGapColCount = 0;
    double simg = 0, simng = 0;
    for(unsigned i = 0; i < uColCount; i++)
    {   
        bool gapColumn = false;
        std::unordered_map <char, unsigned> columnMap;
        for(unsigned j = 0; j < uRowCount; j++)
        {
            if(msa.IsGap(j,i))
            {
                gapColumn = true;
            }
            else
            {
                columnMap[msa.GetChar(j,i)]++;
            }
        }
        auto x = std::max_element(columnMap.begin(), columnMap.end(), [](const std::pair<char, unsigned>& p1, const std::pair<char, unsigned>& p2) {return p1.second < p2.second; });
        unsigned max = x->second;
        if(gapColumn)
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

    for(unsigned i = 0; i < uRowCount ; i++)
    {
        unsigned gapCount = 0;
        unsigned nonGapCount = 0;
        for(unsigned j = 0; j < uColCount; j++)
        {
            if(msa.IsGap(i,j))
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

#endif /* MANUTILITY_H */

