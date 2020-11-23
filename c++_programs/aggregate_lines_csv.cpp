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
    map<string, vector<string>, comparer> aggregated_tissue;
    string header_line;
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
        if (getline(file_in, line)) {
            if (header_line.empty()) header_line = line;
            while (getline(file_in, line)) {
                if (line.length() > 1 and (line[0] == 'S' and line[1] == 'c' or line[0] == 'x' and (line[1] == 'p' or line[1] == 'f'))) {
                    vector<string> parsed = split(line, ',');
                    string pos_id = parsed[0] + "|" + parsed[1] + "|" + parsed[2] + "|" + parsed[3] + "|" + parsed[10];
                    map<string, vector<string> >::iterator it = aggregated_tissue.insert(pair<string, vector<string> > (pos_id, vector<string> ())).first;
                    aggregated_lines.insert(pair<string, string> (pos_id, line));
                    bool found = false;
                    int i = 0;
                    while (i < it->second.size() and not found) {
                        if (it->second[i].compare(current_file_id) == 0) found = true;
                        else i++;
                    }
                    if (not found) it->second.push_back(current_file_id);
                }
            }
        }
        file_in.close();
    }
    map<string, int, comparer> aggregated_totals;
    for (map<string, vector<string> >::iterator it = aggregated_tissue.begin(); it != aggregated_tissue.end(); it++) {
        string aggregated_id = it->second[0];
        for (int i = 1; i < it->second.size(); i++) {
            aggregated_id += ";" + it->second[i];
        }
        map<string, int>::iterator totals_it = aggregated_totals.insert(pair<string, int> (aggregated_id, 0)).first;
        totals_it->second++;
    }
    cout << "TOTAL_unique_pairs," << aggregated_lines.size() << endl;
    for (map<string, int>::iterator it = aggregated_totals.begin(); it != aggregated_totals.end(); it++) cout << it->first << "," << it->second << endl;
    cout << endl;
    cout << header_line << ",amphioxus_tissue" << endl;
    for (map<string, string>::iterator it = aggregated_lines.begin(); it != aggregated_lines.end(); it++) {
        string aggregated_id = aggregated_tissue[it->first][0];
        for (int i = 1; i < aggregated_tissue[it->first].size(); i++) {
            aggregated_id += ";" + aggregated_tissue[it->first][i];
        }
        cout << it->second << "," << aggregated_id << endl;
    }
    
}