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

#include "msa.h"

//using namespace std;

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

SCORE calculateSimgSimngScore(const MSA &msa, double simgWeight,double simngWeight)
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
    simg = simg / gapColCount;
    simng = simng / nonGapColCount;
    SCORE result = simg * simgWeight + simng * simngWeight;
    //printf("c++=%lf, %lf, %lf\n", simg, simng, result);
    return -1 * result;
}


SCORE calculateGapScore(const MSA &msa, double gapWeight)
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
        gapScoreSum += 1.0 * gapCount / nonGapCount;

    }

    return -gapWeight * ( 1.0 * gapScoreSum / uRowCount ) ;

}

#endif /* MANUTILITY_H */

