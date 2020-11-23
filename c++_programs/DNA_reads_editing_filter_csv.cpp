#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <stdlib.h>
using namespace std;


template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}


vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

struct comparer
{
    public:
    bool operator()(const std::string x, const std::string y)
    {
         return x.compare(y) < 0;
    }
};



int main (int argc, char* argv[]) {
    if (argc == 2) {
    string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << "CSV_variants_file genomic_variants_file" << endl;
            return 0;
        }
    }
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << "CSV_variants_file genomic_variants_file" << endl;
        return 1;
    }
    const string CSV_in_path(argv[1]);
    const string DNA_in_path(argv[2]);
    
    ifstream CSV_in;
    ifstream DNA_in;
    DNA_in.open(DNA_in_path.c_str());
    if (not DNA_in.is_open()) {
        cerr << "could not open genomic variants file" << endl;
        return 1;
    }
    map<string, vector<int>, comparer> DNA_vars;
    string DNA_line;
    while (getline(DNA_in, DNA_line)) {
        vector<string> line_parsed = split(DNA_line, '\t');
        string scaf_id = line_parsed[0];
        int pos = atoi(line_parsed[1].c_str());
        map<string, vector<int> >::iterator DNA_it = DNA_vars.insert(pair<string, vector<int> > (scaf_id, vector<int> ())).first;
        DNA_it->second.push_back(pos);
    }
    DNA_in.close();
    CSV_in.open(CSV_in_path.c_str());
    if (not CSV_in.is_open()) {
        cerr << "could not open CSV variants file" << endl;
        return 1;
    }
    string line;
    if (getline(CSV_in, line)) {
        cout << line << endl;
        while (getline(CSV_in, line)) {
            vector<string> line_parsed = split (line, ',');
            string scaf_id = line_parsed[0];
            int pos = atoi(line_parsed[1].c_str());
            map<string, vector<int> >::iterator DNA_it = DNA_vars.find(scaf_id);
            if (DNA_it == DNA_vars.end()) cout << line << endl;
            else {
                bool found = false;
                int i = 0;
                while (not found and i < DNA_it->second.size()) {
                        found = (pos == DNA_it->second[i++]);
                }
                if (not found) cout << line << endl;
            }
        }
    }
    CSV_in.close();
}
