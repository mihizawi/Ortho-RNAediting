#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <stdlib.h>
using namespace std;

typedef struct {
    int CDS_start;
    int CDS_end;
    bool positive_strand;
} gtf_info;

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

void init_gene_table (vector<char> &gene_table) {
    gene_table.resize(64);
    gene_table[0] = 'F';
    gene_table[1] = 'F';
    gene_table[2] = 'L';
    gene_table[3] = 'L';
    gene_table[4] = 'S';
    gene_table[5] = 'S';
    gene_table[6] = 'S';
    gene_table[7] = 'S';
    gene_table[8] = 'Y';
    gene_table[9] = 'Y';
    gene_table[10] = 'X';
    gene_table[11] = 'X';
    gene_table[12] = 'C';
    gene_table[13] = 'C';
    gene_table[14] = 'X';
    gene_table[15] = 'W';
    gene_table[16] = 'L';
    gene_table[17] = 'L';
    gene_table[18] = 'L';
    gene_table[19] = 'L';
    gene_table[20] = 'P';
    gene_table[21] = 'P';
    gene_table[22] = 'P';
    gene_table[23] = 'P';
    gene_table[24] = 'H';
    gene_table[25] = 'H';
    gene_table[26] = 'Q';
    gene_table[27] = 'Q';
    gene_table[28] = 'R';
    gene_table[29] = 'R';
    gene_table[30] = 'R';
    gene_table[31] = 'R';
    gene_table[32] = 'I';
    gene_table[33] = 'I';
    gene_table[34] = 'I';
    gene_table[35] = 'M';
    gene_table[36] = 'T';
    gene_table[37] = 'T';
    gene_table[38] = 'T';
    gene_table[39] = 'T';
    gene_table[40] = 'N';
    gene_table[41] = 'N';
    gene_table[42] = 'K';
    gene_table[43] = 'K';
    gene_table[44] = 'S';
    gene_table[45] = 'S';
    gene_table[46] = 'R';
    gene_table[47] = 'R';
    gene_table[48] = 'V';
    gene_table[49] = 'V';
    gene_table[50] = 'V';
    gene_table[51] = 'V';
    gene_table[52] = 'A';
    gene_table[53] = 'A';
    gene_table[54] = 'A';
    gene_table[55] = 'A';
    gene_table[56] = 'D';
    gene_table[57] = 'D';
    gene_table[58] = 'E';
    gene_table[59] = 'E';
    gene_table[60] = 'G';
    gene_table[61] = 'G';
    gene_table[62] = 'G';
    gene_table[63] = 'G';
}

char translate_codon(char first, char second, char third, const vector<char> &gene_table) {
    int index = 0;
    if (third == 'C' or third == 'c') index += 1*1;
    else if (third == 'A' or third == 'a') index += 1*2;
    else if (third == 'G' or third == 'g') index += 1*3;
    if (second == 'C' or second == 'c') index += 4*1;
    else if (second == 'A' or second == 'a') index += 4*2;
    else if (second == 'G' or second == 'g') index += 4*3;
    if (first == 'C' or first == 'c') index += 16*1;
    else if (first == 'A' or first == 'a') index += 16*2;
    else if (first == 'G' or first == 'g') index += 16*3;
    return gene_table[index];
}

int main (int argc, char* argv[]) {
    vector<char> genetic_code_table;
    init_gene_table(genetic_code_table);
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << "variants_file transcripts_amphi_fasta_file gtf_amphi_file" << endl;
            return 0;
        }
    }
    if (argc != 4) {
            cout << "Usage: " << argv[0] << "variants_file transcripts_amphi_fasta_file gtf_amphi_file" << endl;
            return 1;
    }
    const string var_in_path(argv[1]);
    const string amphi_fasta_in_path(argv[2]);
    const string gtf_in_path(argv[3]);
    ifstream var_in;
    ifstream fasta_in;
    ifstream gtf_in;
    fasta_in.open(amphi_fasta_in_path.c_str());
    if (not fasta_in.is_open()) {
        cerr << "could not open amphiox fasta file" << endl;
        return 1;
    }
    map<string, string, comparer> amphi_trans_seqs;
    string current_trans = "";
    int CDS_start, CDS_end;
    map<string, string>::iterator fas_it;
    string fasta_line;
    while (getline(fasta_in, fasta_line)) {
        if (fasta_line[0] == '>')  {
            if (not current_trans.empty()) {
                fas_it->second.erase(CDS_end, string::npos);
                fas_it->second.erase(0, CDS_start - 1);
            }
            vector<string> fasta_parsed = split(fasta_line, ' ');
            if (fasta_parsed.size() < 3) current_trans = "";
            else {
                current_trans = fasta_parsed[0];
                current_trans.erase(0, 1);
                fas_it = amphi_trans_seqs.insert(pair<string, string > (current_trans, "")).first;
                vector<string> CDS_parsed = split(fasta_parsed[2], '=');
                vector<string> CDS_range_parsed = split(CDS_parsed[1], '-');
                CDS_start = atoi(CDS_range_parsed[0].c_str());
                CDS_end = atoi(CDS_range_parsed[1].c_str());
            }
        }
        else {
            if (not current_trans.empty()) fas_it->second += fasta_line;
        }
    }
    fasta_in.close();
    if (not current_trans.empty()) {
        fas_it->second.erase(CDS_end, string::npos);
        fas_it->second.erase(0, CDS_start - 1);
    }
    fasta_in.close();
    gtf_in.open(gtf_in_path.c_str());
    if (not gtf_in.is_open()) {
        cerr << "could not open human proteome fasta file" << endl;
        return 1;
    }
    map<string, vector<gtf_info>, comparer> amphi_gtf_info;
    string gtf_line;
    while (getline(gtf_in, gtf_line)) {
        vector<string> gtf_parsed = split(gtf_line, '\t');
        if (gtf_parsed[2].compare("CDS") == 0) {
            vector<string> gtf_last_field_parsed = split(gtf_parsed[8], ' ');
            string trans_id = gtf_last_field_parsed[3];
            trans_id.erase(trans_id.length() - 2, 2);
            trans_id.erase(0, 1);
            gtf_info n_info;
            n_info.CDS_start = atoi(gtf_parsed[3].c_str());
            n_info.CDS_end = atoi(gtf_parsed[4].c_str());
            n_info.positive_strand = (gtf_parsed[6][0] == '+');
            map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.insert(pair<string, vector<gtf_info> > (trans_id, vector<gtf_info> ())).first;
            gtf_it->second.push_back(n_info);
        }
    }
    gtf_in.close();
    var_in.open(var_in_path.c_str());
    if (not var_in.is_open()) {
        cerr << "could not open homology variants file" << endl;
        return 1;
    }
    string var_line;
    vector<int> pos_h(3);
    while (getline(var_in, var_line)) {
        if (not var_line.empty()) {
            vector<string> amphi_parsed = split(var_line,'\t');
            vector<string> amphi_trans = split(amphi_parsed[12], ':');
            int scaf_pos = atoi(amphi_parsed[1].c_str());
            for (int t = 0; t < amphi_trans.size(); t++) {
                map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.find(amphi_trans[t]);
                if (gtf_it == amphi_gtf_info.end()) {
                    cerr << "Something went wrong, Unable to find gtf info on transcript " << amphi_trans[t] << endl;
                }
                else {
                    int i = 0;
                    int trans_pos = 0;
                    bool exon_found = false;
                    while (not exon_found and i < gtf_it->second.size()) {
                        if (scaf_pos >= gtf_it->second[i].CDS_start and scaf_pos <= gtf_it->second[i].CDS_end) {
                            exon_found = true;
                            if (gtf_it->second[i].positive_strand) trans_pos += scaf_pos - gtf_it->second[i].CDS_start;
                            else trans_pos += gtf_it->second[i].CDS_end - scaf_pos;
                        }
                        else {
                            trans_pos += gtf_it->second[i].CDS_end - gtf_it->second[i].CDS_start + 1;
                            i++;
                        }
                    }
                    if (exon_found) {
                        map<string, string>::iterator fas_it = amphi_trans_seqs.find(amphi_trans[t]);
                        if (fas_it == amphi_trans_seqs.end()) {
                            cerr << "Something went wrong, Unable to find fasta info on transcript " << amphi_trans[t] << endl; 
                        }
                        else {
                            bool wrong_nucleotide = false;
                            char ori_aa, edi_aa;
                            char first, second, third;
                            if (trans_pos%3 == 0) {
                                first = fas_it->second[trans_pos];
                                second = fas_it->second[trans_pos + 1];
                                third = fas_it->second[trans_pos + 2];
                                if (first != 'A' and first != 'a') {
                                    cerr << "Warning: Something went wrong, the nucleaotide in  " << amphi_trans[t] << " at position " << trans_pos << " is a " << first << ", not an A!" << endl; 
                                    wrong_nucleotide = true;
                                }
                                ori_aa = translate_codon(first, second, third, genetic_code_table);
                                edi_aa = translate_codon('G', second, third, genetic_code_table);
                            }
                            else if (trans_pos%3 == 1) {
                                first = fas_it->second[trans_pos - 1];
                                second = fas_it->second[trans_pos];
                                third = fas_it->second[trans_pos + 1];
                                if (second != 'A' and second != 'a') cerr << "Warning: Something went wrong, the nucleaotide in  " << amphi_trans[t] << " at position " << trans_pos << " is a " << second << ", not an A!" << endl; 
                                ori_aa = translate_codon(first, second, third, genetic_code_table);
                                edi_aa = translate_codon(first, 'G', third, genetic_code_table);
                            }
                            else {
                                first = fas_it->second[trans_pos - 2];
                                second = fas_it->second[trans_pos - 1];
                                third = fas_it->second[trans_pos];
                                if (third != 'A' and third != 'a') cerr << "Warning: Something went wrong, the nucleaotide in  " << amphi_trans[t] << " at position " << trans_pos << " is a " << third << ", not an A!" << endl; 
                                ori_aa = translate_codon(first, second, third, genetic_code_table);
                                edi_aa = translate_codon(first, second, 'G', genetic_code_table);
                            }
                            if (not wrong_nucleotide and ori_aa == edi_aa) cout << var_line << endl;
                        }
                    }
                }
            }
        }
    }
}