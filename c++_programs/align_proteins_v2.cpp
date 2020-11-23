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
    const string HUMAN_PROT_DIR = "human_prot/";
    const string AMPHI_PROT_DIR = "amphi_prot/";
    const string ALIGNMENTS_DIR = "alignments/";
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << "variants_homology_file" << endl;
            return 0;
        }
    }
    if (argc != 2) {
            cout << "Usage: " << argv[0] << "variants_homology_file" << endl;
            return 1;
    }
    const string var_in_path(argv[1]);
    bool debug = false;
    ifstream var_in;
    var_in.open(var_in_path.c_str());
    if (not var_in.is_open()) {
        cerr << "could not open homology variants file" << endl;
        return 1;
    }
    vector<string> var_in_path_parsed = split(var_in_path,'/');
    vector<string> var_in_filename_parsed = split(var_in_path_parsed[var_in_path_parsed.size() - 1],'_');
    string tissue_name = var_in_filename_parsed[2];
    string log_name = tissue_name + "_align.log";
    string all_align_name = tissue_name + "_align_all.fa";
    string log_path = ALIGNMENTS_DIR + tissue_name + "/" + log_name;
    string all_align_path = ALIGNMENTS_DIR + tissue_name + "/" + all_align_name;
    ofstream fall_out(all_align_path.c_str(), ofstream::app);
    if (not fall_out.is_open()) {
        cerr << "could not open file " <<  all_align_path << endl;
        return 1;
    }
    string var_line;
    vector<string> amphi_trans_to_align;;
    vector<string> human_genes_aligned;
    while (getline(var_in, var_line)) {
        if (not var_line.empty()) {
            if (var_line[0] == 'c') {
                string human_id;
                vector<string> human_parsed = split(var_line,'\t');
                if (human_parsed[9].compare("exonic") == 0 and human_parsed[12].compare("UNKNOWN") != 0) { 
                    
                    vector<string> human_field_parsed = split(human_parsed[12], ':');
                    string human_id = human_field_parsed[1];
                    int j = 0;
                    bool gene_found = false;
                    while (j < human_genes_aligned.size() and not gene_found) {
                        if (human_genes_aligned[j++].compare(human_id) == 0) gene_found = true;
                    }
                    if (not gene_found) {
                        human_genes_aligned.push_back(human_id);
                        for (int j = 0; j < amphi_trans_to_align.size(); j++) {
                            ofstream fout(TEMP_FASTA_PATH.c_str(), ofstream::trunc);
                            if (not fout.is_open()) {
                                cerr << "could not open temporary fasta file" <<  TEMP_FASTA_PATH << endl;
                                return 1;
                            }
                            string target_file_path = ALIGNMENTS_DIR + tissue_name + "/";
                            ofstream flog(log_path.c_str(), ofstream::app);
                            if (not flog.is_open()) {
                                cerr << "could not open log file" <<  log_path << endl;
                                return 1;
                            }
                            string amphi_trans_path = AMPHI_PROT_DIR + amphi_trans_to_align[j] + ".fa";
                            ifstream fain(amphi_trans_path.c_str());
                            if (not fain.is_open()) {
                                flog << "WARNING: transcript file " << amphi_trans_path << " was not found" << endl;
                            }
                            else {
                                target_file_path += amphi_trans_to_align[j] + "_";
                                string trans_line;
                                while (getline(fain, trans_line)) fout << trans_line << endl;
                                fain.close();
                            }
                            string human_gene_path = HUMAN_PROT_DIR + human_id + ".fa";
                            ifstream fin(human_gene_path.c_str());
                            if (not fin.is_open()) {
                                flog << "WARNING: transcript file " << human_gene_path << " was not found" << endl;
                                fout.close();
                                flog.close();
                            }
                            else {
                                target_file_path += human_id + "_align.fa";
                                string trans_line;
                                while (getline(fin, trans_line)) fout << trans_line << endl;
                                fin.close();
                                fout.close();
                                flog.close();
                                string command = "mafft " + TEMP_FASTA_PATH + " > " + target_file_path + " 2>> " + log_path;
                                system(command.c_str());
                                ifstream ftarg_in(target_file_path.c_str());
                                if (not ftarg_in.is_open()) {
                                    cerr << "could not open file" << target_file_path << endl;
                                    return 1;
                                }
                                string copy_line;
                                while (getline(ftarg_in, copy_line)) fall_out << copy_line << endl;
                                ftarg_in.close();
                                fall_out << endl;
                            }
                            target_file_path = ALIGNMENTS_DIR + tissue_name + "/";
                        }
                        
                    }
                }
            }
            else {
                vector<string> amphi_parsed = split(var_line,'\t');
                vector<string> amphi_trans = split(amphi_parsed[12], ':');
                for (int t = 0; t < amphi_trans.size(); t++) {
                    int i = 0;
                    bool trans_found = false;;
                    while (i < amphi_trans_to_align.size() and not trans_found) {
                        if (amphi_trans_to_align[i++].compare(amphi_trans[t]) == 0) trans_found = true;
                    }
                    if (not trans_found) amphi_trans_to_align.push_back(amphi_trans[t]);
                }
            }
        }
        else {
            amphi_trans_to_align.clear();
            human_genes_aligned.clear();
        }
    }
    fall_out.close();
    string command_rm = "rm " + TEMP_FASTA_PATH;
    system(command_rm.c_str());
}