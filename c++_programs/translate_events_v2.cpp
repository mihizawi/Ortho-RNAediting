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
    const string TEMP_FASTA_PATH = "temp.fa";
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << "variants_homology_file transcripts_amphi_fasta_file protein_homo_fasta_file" << endl;
            return 0;
        }
    }
    if (argc != 4) {
            cout << "Usage: " << argv[0] << "variants_homology_file transcripts_amphi_fasta_file protein_homo_fasta_file" << endl;
            return 1;
    }
    const string var_in_path(argv[1]);
    const string amphi_fasta_in_path(argv[2]);
    const string human_protein_fasta_in_path(argv[3]);
    bool debug = false;
    ifstream var_in;
    ifstream amphi_fa_in;
    ifstream human_prot_fa_in;
    var_in.open(var_in_path.c_str());
    if (not var_in.is_open()) {
        cerr << "could not open homology variants file" << endl;
        return 1;
    }
    amphi_fa_in.open(amphi_fasta_in_path.c_str());
    if (not amphi_fa_in.is_open()) {
        cerr << "could not open amphiox fasta file" << endl;
        return 1;
    }

    human_prot_fa_in.open(human_protein_fasta_in_path.c_str());
    if (not human_prot_fa_in.is_open()) {
        cerr << "could not open human proteome fasta file" << endl;
        return 1;
    }

    map<string, vector<string>, comparer > amphi_sequences;
    string fasta_line;
    string current_trans = "no trans";
    map<string, vector<string> >::iterator fas_it;
    while (getline(amphi_fa_in, fasta_line)) {
        if (fasta_line[0] == '>') {
            vector<string> fasta_parsed = split(fasta_line, ' ');
            current_trans = fasta_parsed[0];
            current_trans.erase(0, 1);
            fas_it = amphi_sequences.insert(pair<string, vector<string> > (current_trans, vector<string> ())).first;
        }
        fas_it->second.push_back(fasta_line);
    }
    map<string, vector<string>, comparer > human_proteins;
    while (getline(human_prot_fa_in, fasta_line)) {
        if (fasta_line[0] == '>') {
            vector<string> fasta_parsed = split(fasta_line, ' ');
            vector<string> trans_field_parsed = split(fasta_parsed[4], ':');
            current_trans = trans_field_parsed[1];
            fas_it = human_proteins.insert(pair<string, vector<string> > (current_trans, vector<string> ())).first;
        }
        fas_it->second.push_back(fasta_line);
    }
    string var_line;
    while (getline(var_in, var_line)) {
        if (not var_line.empty()) {
            if (var_line[0] == 'c') {
                vector<string> human_parsed = split(var_line,'\t');
                if (human_parsed[9].compare("exonic") == 0 and human_parsed[12].compare("UNKNOWN") != 0) {
                    vector<string> human_field_parsed = split(human_parsed[12], ':');
                    string human_id = human_field_parsed[1];
                    string fasta_file_path = "human_prot/" + human_id + ".fa";
                    ifstream fin(fasta_file_path.c_str());
                    if (fin.good()) fin.close(); 
                    else {
                        ofstream fout(fasta_file_path.c_str());
                        if (not fout.is_open()) {
                            cerr << "could not create file " << fasta_file_path << endl;
                            return 1;
                        }
                        fas_it = human_proteins.find(human_id);
                        if (fas_it != human_proteins.end()) {
                            for (int j = 0; j < fas_it->second.size(); j++) {
                                fout << fas_it->second[j] << endl;
                            }
                        }
                        fout.close();
                    }
                }
            }
            else {
                vector<string> amphi_parsed = split(var_line,'\t');
                vector<string> amphi_trans = split(amphi_parsed[12], ':');
                for (int t = 0; t < amphi_trans.size(); t++) {
                    fas_it = amphi_sequences.find(amphi_trans[t]);
                    if (fas_it != amphi_sequences.end()) {
                        vector<string> first_line_parsed = split(fas_it->second[0], ' ');
                        if (first_line_parsed.size() > 2) {
                            string fasta_file_path = "amphi_prot/" + amphi_trans[t] + ".fa";
                            ifstream fin(fasta_file_path.c_str());
                            if (fin.good()) fin.close(); 
                            else {
                                ofstream fout(TEMP_FASTA_PATH.c_str(), ofstream::trunc);
                                if (not fout.is_open()) {
                                    cerr << "could not open temporary fasta file" <<  TEMP_FASTA_PATH << endl;
                                    return 1;
                                }
                                for (int i = 0; i < fas_it->second.size(); i++) {
                                    fout << fas_it->second[i] << endl;
                                }
                                fout.close();
                                vector<string> CDS_parsed = split(first_line_parsed[2], '=');
                                string command = "transeq " + TEMP_FASTA_PATH + " " + fasta_file_path + " -reg=" + CDS_parsed[1];
                                system(command.c_str());
                            }
                        }
                    }
                }
            }
        }
    }
    string command_rm = "rm " + TEMP_FASTA_PATH;
    system(command_rm.c_str());
}