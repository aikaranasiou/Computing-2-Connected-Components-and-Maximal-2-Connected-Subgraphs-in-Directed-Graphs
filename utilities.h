/*
 Copyright (C) 2015 Nilakantha Paudel, Loukas Georgiadis, Giuseppe F. Italiano
 All rights reserved.

 This software may be modified and distributed under the terms
 of the BSD license as followoing.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation  and/or other materials provided with the distribution.

 3. Neither the names of the copyright holders nor the names of any
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
 */

//============================================================================
// Project     : TWO VERTEX STRONGLY CONNECTED COMPONENT
// Name		   : utilities.h
// Author      : Nilakantha Paudel [Nilu]
// Version     : 1.0
// Description : It contains the general utitlity methods for the project.
//============================================================================
#ifndef TWOVSCC_COMMON_UTILITIES_H_
#define TWOVSCC_COMMON_UTILITIES_H_

#include "common_system_header.h"

class Utilities final {

public:
    //public Typedefs and Enums -- not defined.
    //public static const data members -- not defined.
    static const int LOG_BASE = 2;

    //public Constructors and assignment operators

    Utilities() {
    }
    //Destructor
    virtual ~Utilities() {
    }

    template<typename T>
    static string TO_STRING(T val) {
        stringstream stream;
        stream << val;
        return stream.str();
    }

    //public Methods, including static methods
    static int GetLOGwithBase(int param, int base) {
        return ceil(((double) log2(param) / (double) log2(base)));
    }
    static int GetPOWERwithBase(int base, int power) {
        return (int) (pow(base, power));
    }
    static string GetCurrentDateTime() {

        time_t now = time(0);
        struct tm tstruct;
        char buf[80];
        tstruct = *localtime(&now);
        strftime(buf, sizeof(buf), "%Y/%m/%d @ %X", &tstruct);
        return buf;

    }
    static string GetFileName(string FileName) {
        string output_file_name;
        vector<string> vList = SplitString(".", FileName);
        int size = (int) vList.size();
        for (int i = 0; i < (size - 1); i++) {
            output_file_name += vList.at(i);
        }
        return output_file_name;
    }
    // trim string from both ends
    static string GetRatioString(int n, int m) {

        float val = (float) (m / (float) n);
        string finalString = Utilities::TO_STRING(val);
        return Utilities::TrimString(finalString);
    }
    static std::string &TrimString(std::string &s) {
        return Utilities::LTRIM(Utilities::RTRIM(s));
    }
    static vector<string> GetSplitedStringByDelimeter(string inputString,
            string delimiter) {

        vector<string> finalList;

        size_t pos = 0;
        std::string token;
        while ((pos = inputString.find(delimiter)) != std::string::npos) {
            token = inputString.substr(0, pos);
            //std::cout << token << std::endl;
            finalList.push_back(token);
            inputString.erase(0, pos + delimiter.length());
        }
        //std::cout << inputString << std::endl;
        finalList.push_back(inputString);

        return finalList;
    }

    static void SortTuples(int row_size, int column_size, int *tuples) {

        int total_elements = row_size * column_size;

        int *temp_tuple = new int[total_elements];
        int* rows = new int[row_size + 2];

        for (int column_index = (column_size - 1); column_index >= 0;
                column_index--) {

            for (int j = 0; j <= row_size; j++)
                rows[j] = 0; //initialize

            //count how many times every elements apears in that column (columnIndex).
            for (int e = 0; e < row_size; e++)
                rows[tuples[column_size * e + column_index] + 1]++;

            // cumulative value of the counts for that column (columnIndex)
            for (int j = 1; j <= row_size; j++)
                rows[j] += rows[j - 1];

            //sort the row according to the elements of that column (columnIndex)
            // And then put int right array
            for (int row_index = 0; row_index < row_size; row_index++) {

                int pre_row_index = row_index * column_size + column_index;

                int k = rows[tuples[pre_row_index]];

                int new_row_index = k * column_size;

                for (int e = 0; e < column_size; e++)
                    temp_tuple[new_row_index + e] = tuples[column_size
                            * row_index + e];

                // increase the index for the row of that element
                rows[tuples[pre_row_index]]++;
            }

            //copy the triples to the original array
            for (int j = 0; j < total_elements; j++)
                tuples[j] = temp_tuple[j];
        }
        delete[] temp_tuple;
        delete[] rows;
    }

    //Public Data Members

private:
    //private Typedefs and Enums -- not defined.
    //private static const data members -- not defined.

    //private Methods, including static methods

    static vector<string> SplitString(const string &delimiter,
            const string &inputString) {

        vector<string> arr;

        int strleng = inputString.length();
        int delleng = delimiter.length();
        if (delleng == 0)
            return arr; //no change

        int i = 0;
        int k = 0;
        while (i < strleng) {
            int j = 0;
            while (i + j < strleng && j < delleng
                    && inputString[i + j] == delimiter[j])
                j++;
            if (j == delleng) //found delimiter
                    {
                arr.push_back(inputString.substr(k, i - k));
                i += delleng;
                k = i;
            } else {
                i++;
            }
        }
        arr.push_back(inputString.substr(k, i - k));
        return arr;
    }
    // trim from start
    static std::string &LTRIM(std::string &s) {
        s.erase(s.begin(),
                std::find_if(s.begin(), s.end(),
                        std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    // trim from end
    static std::string &RTRIM(std::string &s) {
        s.erase(
                std::find_if(s.rbegin(), s.rend(),
                        std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
                s.end());
        return s;
    }

    //private Data Members

};
#endif /* TWOVSCC_COMMON_UTILITIES_H_ */
