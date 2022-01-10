// SignificantFigures.c
// David Grund, Oct 27, 2021

#include <stdio.h> // printf
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()

// https://thispointer.com/c-convert-double-to-string-and-manage-precision-scientific-notation/
// https://www.quora.com/How-do-I-determine-the-number-of-digits-after-decimal-in-a-floating-point-number-in-C++

Int_t NoDecPlaces(Double_t d);

void SignificantFigures(){

    Double_t d = 10504.2548;

    ofstream outfile("out.txt");
    outfile << std::fixed << std::setprecision(1);
    outfile << d;
    outfile << "\n";
    outfile << std::defaultfloat << std::setprecision(1);
    outfile << d;
    outfile.close();

    Printf("%i", NoDecPlaces(41548.000045451));

    return;
}

///*
Int_t NoDecPlaces(Double_t d){
    stringstream ss;
    ss << std::fixed << std::setprecision(10) << d;
    string s;
    ss >> s;
    Printf("%s", s.data());
    if(s[0] != "0"){
        Int_t i = 0;
        while(s[i] != ".") i++;
        if(i == 1)
    }
    else 
    if(s.length() == 1) return 1;
    else return s.length() - 2;
    return 1;
}
//*/