#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
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
            cout << "Usage: " << argv[0] << "RNAediting_homo blast_amphi_homo" << endl;
            return 0;
        }
    }
    if (argc != 3) {
            cout << "Usage: " << argv[0] << "RNAediting_homo blast_amphi_homo" << endl;
            return 1;
    }
    const string human_in_path(argv[1]);
    const string blast_in_path(argv[2]);
    bool debug = false;
    ifstream human_in;
    ifstream blast_in;
    blast_in.open(blast_in_path.c_str());
    if (not blast_in.is_open()) {
        cerr << "could not open blast results file" << endl;
        return 1;
    }
    human_in.open(human_in_path.c_str());
    if (not human_in.is_open()) {
        cerr << "could not open human editing file" << endl;
        return 1;
    }
    vector<string> human_editing;
    string human_line;
    while (getline(human_in, human_line)) {
        human_editing.push_back(human_line);
    }
    string blast_line;
    string current_gene = "no_gene";
    vector<string> homologs_found;
    map<string, vector<vector<int> >, comparer> homolog_editing;
    while (getline(blast_in, blast_line)) {
        if (blast_line[0] != '#') {
            vector<string> blast_parsed = split(blast_line, '\t');
            vector<string> trans_id_parsed = split(blast_parsed[0], '_');
            string amphi_gene = trans_id_parsed[0];
            if (current_gene.compare(amphi_gene) != 0) {
                map<string, vector< vector<int> > >::iterator it;
                it = homolog_editing.find(current_gene);
                if (it != homolog_editing.end()) {
                    int size1 = it->second.size();
                    cout << current_gene << " " << size1;
                    for (int i = 0; i < size1; i++) {
                        int size2 = it->second[i].size();
                        cout << " " << size2;
                        for (int j = 0; j < size2; j++) {
                            cout << " " << it->second[i][j];
                        }
                    }
                    cout << endl;
                }
                homologs_found.clear();
                current_gene = amphi_gene;
            }
            if (homologs_found.size() < 5) {
                vector<string> human_id_parsed = split(blast_parsed[1], '|');
                string human_gene = human_id_parsed[2];
                bool seen_homolog = false;
                for (int i = 0; i < homologs_found.size(); i++) {
                    if (homologs_found[i].compare(human_gene) == 0) seen_homolog = true;
                }
                if (not seen_homolog) {
                    homologs_found.push_back(human_gene);
                    map<string, vector< vector<int> > >::iterator it;
                    it = homolog_editing.find(amphi_gene);
                    if (it == homolog_editing.end()) {
                        vector<vector<int> > v;
                        v.push_back(vector<int> ());
                        for (int i = 0; i < human_editing.size(); i++) {
                            vector<string> human_parsed = split(human_editing[i],'\t');
                            if (human_parsed[9].compare("exonic") == 0 and human_parsed[10].compare(human_gene) == 0) {
                                v[0].push_back(i);
                            }
                        }
                        if (v[0].size() > 0) {
                            homolog_editing.insert(pair<string, vector< vector<int> > > (amphi_gene, v));
                        }
                    }
                    else {
                        vector<int> v;
                        for (int i = 0; i < human_editing.size(); i++) {
                            vector<string> human_parsed = split (human_editing[i],'\t');
                            if (human_parsed[9].compare("exonic") == 0 and human_parsed[10].compare(human_gene) == 0) {
                                v.push_back(i);
                            }
                        }
                        
                        if (v.size() > 0) {
                            it->second.push_back(v);
                        }
                    }
                }
            }
        }
    }
}