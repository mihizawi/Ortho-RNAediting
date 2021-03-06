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
            cout << "Usage: " << argv[0] << "[options] [input files]" << endl;
            cout << "Options: -p file prefix -s file suffx" << endl;
            return 0;
        }
    }
    string file_prefix;
    string file_suffix;
    vector<string> input_paths;
    for (int i = 1; i < argc; i++) {
        if (string(argv[i]) == "-p") {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << "[options] [input files]" << endl;
                cerr << "Options: -p file prefix -s file suffx" << endl;
                return 1;
            }
            file_prefix = argv[i];
        }
        else if (string(argv[i]) == "-s") {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << "[options] [input files]" << endl;
                cerr << "Options: -p file prefix -s file suffx" << endl;
                return 1;
            }
            file_suffix = argv[i];
        }
        else input_paths.push_back(argv[i]);
    }
    if (input_paths.empty()) {
        cerr << "Usage: " << argv[0] << "[options] [input files]" << endl;
        cerr << "Options: -p file prefix -s file suffx" << endl;
        return 1;
    }
    map<string, string, comparer> aggregated_lines;
    for (int i = 0; i < input_paths.size(); i++) {
        ifstream file_in;
        file_in.open(input_paths[i].c_str());
        if (not file_in.is_open()) {
            cerr << "could not open input file " << input_paths[i] << endl;
            return 1;
        }
        string current_file_id;
        if (input_paths[i].find('/') == string::npos) current_file_id = input_paths[i];
        else {
            vector<string> path_split = split(input_paths[i], '/');
            current_file_id = path_split[path_split.size() - 1];
        }
        size_t pos_prefix = current_file_id.find(file_prefix);
        if (pos_prefix != string::npos) current_file_id.erase(0, pos_prefix + file_prefix.length());
        size_t pos_suffix = current_file_id.find(file_suffix);
        if (pos_suffix != string::npos) current_file_id.erase(current_file_id.begin() + pos_suffix, current_file_id.end());
        string line;
        while (getline(file_in, line)) {
            if (line.length() > 1 and (line[0] == 'S' and line[1] == 'c' or line[0] == 'x' and (line[1] == 'p' or line[1] == 'f'))) {
                vector<string> parsed = split(line, '\t');
                string pos_id = parsed[0] + "|" + parsed[1];
                aggregated_lines.insert(pair<string, string> (pos_id, line)).first;
            }
        }
        file_in.close();
    }
    for (map<string, string>::iterator it = aggregated_lines.begin(); it != aggregated_lines.end(); it++) {
        cout << it->second << endl;
    }
    
}