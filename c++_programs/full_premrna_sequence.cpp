#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <stdlib.h>
#include <map>
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

typedef struct {
    string gene_id;
    string trans_id;
    int gtf_start;
    int gtf_end;
    bool positive_strand;
    bool CDS;
} gtf_info;

struct comparer
{
    public:
    bool operator()(const std::string x, const std::string y)
    {
         return x.compare(y) < 0;
    }
};

string complementary_seq(const string &seq) {
    string cmp(seq);
    char c;
    for (int i = 0; i < seq.length(); i++) {
        c = seq[i];
        if (c == 'A') cmp[i] = 'T';
        else if (c == 'a') cmp[i] = 't';
        else if (c == 'T') cmp[i] = 'A';
        else if (c == 't') cmp[i] = 'a';
        else if (c == 'G') cmp[i] = 'C';
        else if (c == 'g') cmp[i] = 'c';
        else if (c == 'C') cmp[i] = 'G';
        else if (c == 'c') cmp[i] = 'g';
    }
    return cmp;
}

int main (int argc, char* argv[]) {
    if (argc == 2) {
    string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << "filtered_variants_size genome_fasta_file gtf_file [options]" << endl;
            cout << "Usage: " << " options: -c (ignore noncoding variants)      -i (ignore nongenic variants)" << endl;
            return 0;
        }
    }
    if (argc < 3) {
            cerr << "Usage: " << argv[0] << "filtered_variants_size genome_fasta_file gtf_file [options]" << endl;
            cerr << "Usage: " << " options: -c (ignore noncoding variants)      -i (ignore nongenic variants)" << endl;
            return 1;
    }
    const string var_in_path(argv[1]);
    const string fasta_in_path(argv[2]);
    const string gtf_in_path(argv[3]);
    bool ignore_nongene = false;
    bool ignore_noncoding = false;
    for (int i = 4; i < argc; i++) {
        if (string(argv[i]) == "-c") {
            ignore_noncoding = true;
        }
        else if (string(argv[i]) == "-i") {
            ignore_nongene = true;
        }        
        else {
            cerr << "Usage: " << argv[0] << "filtered_variants_size genome_fasta_file gtf_file [options]" << endl;
            cerr << "Usage: " << " options: -c (ignore noncoding variants)      -i (ignore nongenic variants)" << endl;
            return 1;
        }
    }
    ifstream var_in;
    ifstream fasta_in;
    var_in.open(var_in_path.c_str());
    if (not var_in.is_open()) {
        cerr << "could not open filtered variants file" << endl;
        return 1;
    }
    fasta_in.open(fasta_in_path.c_str());
    if (not fasta_in.is_open()) {
        cerr << "could not open fasta file" << endl;
        return 1;
    }
    string fasta_line; 
    getline(fasta_in, fasta_line);
    if (fasta_line[0] != '>') {
        cerr << "first line in fasta file is not a valid scaffold id" << endl;
        return 1;
    }
    fasta_line.erase(0, 1);
    vector<string> fasta_parsed = split(fasta_line, ' ');
    string scaf_id = fasta_parsed[0];
    string prev_scaf_id = fasta_parsed[0];
    fasta_parsed[1].erase(0, 5);
    vector<string> size_parsed = split(fasta_parsed[1], 'b');
    int scaf_size = atoi(size_parsed[0].c_str());
    string scaf_seq = string();
    scaf_seq.reserve(scaf_size);
    bool eof = false;
    bool end_scaffold = false;
    while (not end_scaffold and not eof) {
        if (getline(fasta_in, fasta_line)) {
            if (fasta_line[0] != '>') {
                scaf_seq += fasta_line;
            }
            else {
                end_scaffold = true;
                fasta_line.erase(0, 1);
                fasta_parsed = split(fasta_line, ' ');
                scaf_id = fasta_parsed[0];
            }
        }
        else {
            eof = true;
        }
    }
    ifstream gtf_in;
    gtf_in.open(gtf_in_path.c_str());
    if (not gtf_in.is_open()) {
        cerr << "could not open transcript reference file" << endl;
        return 1;
    }
    map<string, vector<gtf_info>, comparer> amphi_gtf_info;
    string gtf_line;
    while (getline(gtf_in, gtf_line)) {
        vector<string> gtf_parsed = split(gtf_line, '\t');
        bool is_coding = (gtf_parsed[2].compare("CDS") == 0);
        if (is_coding or not (ignore_nongene and ignore_noncoding)) {
            string scaf_id = gtf_parsed[0];
            gtf_info n_info;
            vector<string> gtf_last_field_parsed = split(gtf_parsed[8], ' ');
            n_info.gene_id = gtf_last_field_parsed[1];
            n_info.gene_id.erase(n_info.gene_id.length() - 2, 2);
            n_info.gene_id.erase(0, 1);
            n_info.trans_id = gtf_last_field_parsed[3];
            n_info.trans_id.erase(n_info.trans_id.length() - 2, 2);
            n_info.trans_id.erase(0, 1);
            n_info.gtf_start = atoi(gtf_parsed[3].c_str());
            n_info.gtf_end = atoi(gtf_parsed[4].c_str());
            n_info.positive_strand = (gtf_parsed[6][0] == '+');
            n_info.CDS = is_coding;
            map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.insert(pair<string, vector<gtf_info> > (scaf_id, vector<gtf_info> ())).first;
            gtf_it->second.push_back(n_info);
        }

    }
    string line;
    while (getline(var_in, line)) {
        vector<string> parsed = split (line, '\t');
        if ((not ignore_nongene or parsed[11].compare("unknown_gene") != 0) and (not ignore_noncoding or parsed[12].compare("coding") == 0)) { 
            if (parsed[0].compare(prev_scaf_id) != 0) {
                while (parsed[0].compare(scaf_id) != 0 and not eof) {
                    if (getline(fasta_in, fasta_line)) {
                        if (fasta_line[0] == '>') {
                            fasta_line.erase(0, 1);
                            fasta_parsed = split(fasta_line, ' ');
                            scaf_id = fasta_parsed[0];
                        }
                    }
                    else eof = true;
                }
                if (eof) {
                    fasta_in.close();
                    fasta_in.open(fasta_in_path.c_str());
                    if (not fasta_in.is_open()) {
                        cerr << "could not open fasta file" << endl;
                        return 1;
                    }
                    getline(fasta_in, fasta_line);
                    fasta_line.erase(0, 1);
                    fasta_parsed = split(fasta_line, ' ');
                    scaf_id = fasta_parsed[0];
                    eof = false;
                }
                while (parsed[0].compare(scaf_id) != 0 and not eof) {
                    if (getline(fasta_in, fasta_line)) {
                        if (fasta_line[0] == '>') {
                            fasta_line.erase(0, 1);
                            fasta_parsed = split(fasta_line, ' ');
                            scaf_id = fasta_parsed[0];
                        }
                    }
                    else eof = true;
                }
                if (eof) {
                    cerr << "scaffold " << parsed[0] << " does not match smy scaffold in the fasta file" << endl;
                    return 1;
                }
                prev_scaf_id = fasta_parsed[0];
                scaf_seq.clear();
                if (line[0] == 'x') {
                    vector<string> field_parsed = split (fasta_parsed[1], ';');
                    vector<string> length_parsed = split (field_parsed[3], ':');
                    scaf_size = atoi(length_parsed[1].c_str());
                }
                else {
                    fasta_parsed[1].erase(0, 5);
                    size_parsed = split(fasta_parsed[1], 'b');
                    scaf_size = atoi(size_parsed[0].c_str());
                }
                scaf_seq.reserve(scaf_size + scaf_size/5);
                end_scaffold = false;
                while (not end_scaffold and not eof) {
                    if (getline(fasta_in, fasta_line)) {
                        if (fasta_line[0] != '>') {
                            scaf_seq += fasta_line;
                        }
                        else {
                            end_scaffold = true;
                            fasta_line.erase(0, 1);
                            fasta_parsed = split(fasta_line, ' ');
                            scaf_id = fasta_parsed[0];
                        }
                    }
                    else {
                        eof = true;
                    }
                }
            }
            scaf_size = scaf_seq.length();
            map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.find(parsed[0]);
            if (gtf_it == amphi_gtf_info.end()) {
                cerr << "Something went wrong, Unable to find scaffold " << parsed[0] << " in gtf file" << endl; 
                return 1;
            }
            vector<string> trans_ids = split(parsed[12], ':');
            for (int i = 0; i < trans_ids.size(); i++) {
                string current_t_id = trans_ids[i];
                int start = scaf_size;
                int end = 0;
                for (int j = 0; j < gtf_it->second.size(); j++) {
                    if (gtf_it->second[j].trans_id.compare(current_t_id) == 0) {
                        if (gtf_it->second[j].gtf_start < start) start = gtf_it->second[j].gtf_start;
                        if (gtf_it->second[j].gtf_end > end) end = gtf_it->second[j].gtf_end;
                    }
                }
                int pos = atoi(parsed[1].c_str());
                int gene_length = end - start + 1;
                string seq(scaf_seq, start - 1, gene_length);
                if (parsed[10][0] == '-') seq = complementary_seq(seq);
                int trans_pos = pos - start + 1;
                cout << ">" << current_t_id << '\t' << trans_pos << '\t' << line << endl;
                //cout << start << " " << end << " " << gene_length << endl;
                cout << seq << endl;
            }
        }
    }
    var_in.close();
    fasta_in.close();
}