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

typedef struct {
    string gene_id;
    string trans_id;
    int position;
    int aligned_position;
    char ori_aa;
    char edi_aa;
    string line;
} editing_info;

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
            cout << "Usage: " << argv[0] << " variants_homology_file amphi_transcriptome_fasta_file amphi_gtf_file [options]" << endl;
            cout << "Options: -r (aminoacid position difference range) -a (alignments directory path) -p (output directory path)" << endl;
            return 0;
        }
    }
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " variants_homology_file amphi_transcriptome_fasta_file amphi_gtf_file [options]" << endl;
        cerr << "Options: -r (aminoacid position difference range) -a (alignments directory path) -p (output directory path)" << endl;
        return 1;
    }
    string ALIGNMENTS_DIR = "alignments/";
    string CLASSIFIED_EDITING_DIR = "classified_editing/";
    const string var_in_path(argv[1]);
    const string fasta_in_path(argv[2]);
    const string gtf_in_path(argv[3]);
    int aaRange = 0;
    for (int i = 4; i < argc; i++) {
        if (string(argv[i]) == "-r") {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << " variants_homology_file amphi_transcriptome_fasta_file amphi_gtf_file [options]" << endl;
                cerr << "Options: -r (aminoacid position difference range) -a (alignments directory path) -p (output directory path)" << endl;
                return 1;
            }
            aaRange = atoi(argv[i]);
        }
        else if (string(argv[i]) == "-a") {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << " variants_homology_file amphi_transcriptome_fasta_file amphi_gtf_file [options]" << endl;
                cerr << "Options: -r (aminoacid position difference range) -a (alignments directory path) -p (output directory path)" << endl;
                return 1;
            }
            const string aldir(argv[i]);
            ALIGNMENTS_DIR = aldir + "/";
        }
        else if (string(argv[i]) == "-p") {
            i++;
            if (i >= argc) {
                cerr << "Usage: " << argv[0] << " variants_homology_file amphi_transcriptome_fasta_file amphi_gtf_file [options]" << endl;
                cerr << "Options: -r (aminoacid position difference range) -a (alignments directory path) -p (output directory path)" << endl;
                return 1;
            }
            const string outPath(argv[i]);
            CLASSIFIED_EDITING_DIR = outPath + "/";
        }
        else {
            cerr << "Usage: " << argv[0] << " variants_homology_file amphi_transcriptome_fasta_file amphi_gtf_file [options]" << endl;
            cerr << "Options: -r (aminoacid position difference range) -a (alignments directory path) -p (output directory path)" << endl;
            return 1;
        }
    }
	bool debug = false;
    ifstream var_in;
    var_in.open(var_in_path.c_str());
    if (not var_in.is_open()) {
        cerr << "could not open homology variants file" << endl;
        return 1;
    }
    vector<string> var_in_path_parsed = split(var_in_path,'/');
    vector<string> var_in_filename_parsed = split(var_in_path_parsed[var_in_path_parsed.size() - 1],'_');
    const string tissue_name = var_in_filename_parsed[2];
    const string tissue_align_dir = ALIGNMENTS_DIR + tissue_name + "/";
    const string tissue_class_dir = CLASSIFIED_EDITING_DIR + tissue_name + "/";
    ifstream fasta_in;
    fasta_in.open(fasta_in_path.c_str());
    if (not fasta_in.is_open()) {
        cerr << "could not open amphi transcript fasta file" << endl;
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
    ifstream gtf_in;
    gtf_in.open(gtf_in_path.c_str());
    if (not gtf_in.is_open()) {
        cerr << "could not open amphi transcript gtf file" << endl;
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
    string var_line;
    vector<editing_info> human_editing;
    vector<editing_info> amphi_editing;
    vector<string> amphi_trans_ids;
    vector<string> human_genes;
    const string out_class_0_path = tissue_class_dir + "editing_homology_" + tissue_name + "_class_0_syn";
    ofstream class_0_out(out_class_0_path.c_str());
    if (not class_0_out.is_open()) {
        cerr << "could not open class 0 output file" << out_class_0_path << endl;
        return 1;
    }
    const string out_class_1_path = tissue_class_dir + "editing_homology_" + tissue_name + "_class_1_EtoE";
    ofstream class_1_out(out_class_1_path.c_str());
    if (not class_1_out.is_open()) {
        cerr << "could not open class 1 output file" << out_class_1_path << endl;
        return 1;
    }
    const string out_class_2_path = tissue_class_dir + "editing_homology_" + tissue_name + "_class_2_NEtoE";
    ofstream class_2_out(out_class_2_path.c_str());
    if (not class_2_out.is_open()) {
        cerr << "could not open class 2 output file" << out_class_2_path << endl;
        return 1;
    }
    const string out_class_3_path = tissue_class_dir + "editing_homology_" + tissue_name + "_class_3_EtoNE";
    ofstream class_3_out(out_class_3_path.c_str());
    if (not class_3_out.is_open()) {
        cerr << "could not open class 3 output file" << out_class_3_path << endl;
        return 1;
    }
    const string out_class_4_path = tissue_class_dir + "editing_homology_" + tissue_name + "_class_4_NEtoNE";
    ofstream class_4_out(out_class_4_path.c_str());
    if (not class_4_out.is_open()) {
        cerr << "could not open class 4 output file" << out_class_4_path << endl;
        return 1;
    }
/*	const string out_class_5_path = tissue_class_dir + "editing_homology_" + tissue_name + "_class_5_DifPos";
	ofstream class_5_out;
	if (aaRange > 0) {
		class_5_out.open(out_class_5_path);
        if (not class_5_out.is_open()) {
            cerr << "could not open class 4 output file" << out_class_5_path << endl;
            return 1;
        }
	}*/
    const string out_aligned_path = tissue_class_dir + "editing_homology_" + tissue_name + "_aligned";
    ofstream aligned_out(out_aligned_path.c_str());
    if (not aligned_out.is_open()) {
        cerr << "could not open aligned output file" << out_aligned_path << endl;
        return 1;
    }
    while (getline(var_in, var_line)) {
        if (var_line.empty()) {
            cerr << endl;
            if (human_editing.size() > 0) {
                for (int i = 0; i < amphi_trans_ids.size(); i++) {
                    string current_amphi_id = amphi_trans_ids[i];
                    for (int g = 0; g < human_genes.size(); g++) {
                        string human_trans_id = human_genes[g];
                        string align_in_path = tissue_align_dir + current_amphi_id + "_" + human_trans_id + "_align.fa";
                        ifstream align_in(align_in_path.c_str());
                        if (not align_in.is_open()) {
                            cerr << "Warning: could not open alignment file " << align_in_path << endl;
                        }
                        else if (align_in.peek() == ifstream::traits_type::eof()) {
                            cerr << "Warning: alignment file " << align_in_path << " is empty" << endl;
                        }
                        else {
                            string align_line;
                            vector<string> aligned_fasta_lines;
                            vector<string> order_of_ids;
                            while (getline(align_in, align_line)) {
                                if (align_line[0] == '>') {
                                    vector<string> align_parsed = split(align_line, ' ');
                                    string current_id;
                                    if (align_line[1] == 'B' and align_line[2] == 'L') {
                                        current_id = align_parsed[0];
                                        current_id.erase(0, 1);
                                        current_id.erase(current_id.length() - 2, 2);
                                    }
                                    else {
                                        vector<string> transcript_parsed = split(align_parsed[4], ':');
                                        current_id = transcript_parsed[1];
                                    }
                                    order_of_ids.push_back(current_id);
                                    aligned_fasta_lines.push_back(align_line);
                                    aligned_fasta_lines.push_back("");
                                }
                                else {
                                    aligned_fasta_lines[aligned_fasta_lines.size() - 1] += align_line;
                                }
                            }
                            align_in.close();
                            int amphi_order = 0;
                            bool id_found = false;
                            while (not id_found and amphi_order < order_of_ids.size()) {
                                if (order_of_ids[amphi_order].compare(current_amphi_id) == 0) id_found = true;
                                else amphi_order++;
                            }
                            for (int j = 0; j < amphi_editing.size(); j++) {
                                if (amphi_editing[j].trans_id.compare(current_amphi_id) == 0) {
                                    int k = 0;
                                    int count = 0;
                                    while (count < amphi_editing[j].position) {
                                        if (aligned_fasta_lines[amphi_order*2 + 1][k++] != '-') count++;
                                    }
                                    amphi_editing[j].aligned_position = k;
                                }
                            }
                            for (int j = 0; j < human_editing.size(); j++) {
                                if (human_editing[j].trans_id.compare(human_trans_id) == 0) {
                                    int order = 0;
                                    bool id_found = false;
                                    while (not id_found and order < order_of_ids.size()) {
                                        if (order_of_ids[order].compare(human_editing[j].trans_id) == 0) id_found = true;
                                        else order++;
                                    }
                                    if (not id_found) {
                                        cerr << "Warning: The protein with transcript id " << human_editing[j].trans_id << " was not found in the alignment file" << endl;
                                        human_editing.erase(human_editing.begin() + j);
                                        j--;
                                    }
                                    else {
                                        int k = 0;
                                        int count = 0;
                                        while (count < human_editing[j].position) {
                                            if (aligned_fasta_lines[order*2 + 1][k++] != '-') count++;
                                        }
                                        human_editing[j].aligned_position = k;
                                    }
                                }
                            }
                            vector<string> modified_alignment_lines(aligned_fasta_lines);
                            for (int j = 0; j < amphi_editing.size(); j++) {
                                if (amphi_editing[j].trans_id.compare(current_amphi_id) == 0) {
                                    vector<int> hjs;
                                    for (int hj = 0; hj < human_editing.size(); hj++) {
                                        if (human_editing[hj].trans_id.compare(human_trans_id) == 0) {
                                            if (amphi_editing[j].aligned_position - human_editing[hj].aligned_position >= -aaRange
                                                and amphi_editing[j].aligned_position - human_editing[hj].aligned_position <= aaRange) hjs.push_back(hj);
                                        }
                                    }
                                    
                                    for (int hjj = 0; hjj < hjs.size(); hjj++) {
                                        /*if (amphi_editing[j].aligned_position != human_editing[hjs[hjj]].aligned_position) {
                                            class_5_out << amphi_editing[j].line << '\t' << amphi_editing[j].ori_aa << amphi_editing[j].position << amphi_editing[j].edi_aa << endl;
                                            class_5_out << human_editing[hjs[hjj]].line << endl;
                                            if (amphi_editing[j].aligned_position - 1 == aligned_fasta_lines[amphi_order*2 + 1].length()) modified_alignment_lines[amphi_order*2 + 1].insert(aligned_fasta_lines[amphi_order*2 + 1].length(), "*5*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < amphi_editing[j].position) {
                                                    if (modified_alignment_lines[amphi_order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[amphi_order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[amphi_order*2 + 1].insert(k, "*5*");
                                            }
                                            int order = 0;
                                            bool id_found = false;
                                            while (not id_found and order < order_of_ids.size()) {
                                                if (order_of_ids[order].compare(human_editing[hjs[hjj]].trans_id) == 0) id_found = true;
                                                else order++;
                                            }
                                            if (human_editing[hjs[hjj]].aligned_position - 1 == aligned_fasta_lines[order*2 + 1].length()) modified_alignment_lines[order*2 + 1].insert(aligned_fasta_lines[order*2 + 1].length(), "*5*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < human_editing[hjs[hjj]].position) {
                                                    if (modified_alignment_lines[order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[order*2 + 1].insert(k, "*5*");
                                            }
                                        }
                                        else*/ if (amphi_editing[j].ori_aa == amphi_editing[j].edi_aa and human_editing[hjs[hjj]].ori_aa == human_editing[hjs[hjj]].edi_aa 
                                             and amphi_editing[j].ori_aa == human_editing[hjs[hjj]].ori_aa) {
                                            class_0_out << amphi_editing[j].line << '\t' << amphi_editing[j].ori_aa << amphi_editing[j].position << amphi_editing[j].edi_aa << endl;
                                            class_0_out << human_editing[hjs[hjj]].line << endl;
                                            if (amphi_editing[j].aligned_position - 1 == aligned_fasta_lines[amphi_order*2 + 1].length()) modified_alignment_lines[amphi_order*2 + 1].insert(aligned_fasta_lines[amphi_order*2 + 1].length(), "*0*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < amphi_editing[j].position) {
                                                    if (modified_alignment_lines[amphi_order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[amphi_order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[amphi_order*2 + 1].insert(k, "*0*");
                                            }
                                            int order = 0;
                                            bool id_found = false;
                                            while (not id_found and order < order_of_ids.size()) {
                                                if (order_of_ids[order].compare(human_editing[hjs[hjj]].trans_id) == 0) id_found = true;
                                                else order++;
                                            }
                                            if (human_editing[hjs[hjj]].aligned_position - 1 == aligned_fasta_lines[order*2 + 1].length()) modified_alignment_lines[order*2 + 1].insert(aligned_fasta_lines[order*2 + 1].length(), "*0*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < human_editing[hjs[hjj]].position) {
                                                    if (modified_alignment_lines[order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[order*2 + 1].insert(k, "*0*");
                                            }
                                        }
                                        else if (amphi_editing[j].ori_aa == human_editing[hjs[hjj]].ori_aa and amphi_editing[j].edi_aa == human_editing[hjs[hjj]].edi_aa) {
                                            class_1_out << amphi_editing[j].line << '\t' << amphi_editing[j].ori_aa << amphi_editing[j].position << amphi_editing[j].edi_aa << endl;
                                            class_1_out << human_editing[hjs[hjj]].line << endl;
                                            if (amphi_editing[j].aligned_position - 1 == aligned_fasta_lines[amphi_order*2 + 1].length()) modified_alignment_lines[amphi_order*2 + 1].insert(aligned_fasta_lines[amphi_order*2 + 1].length(), "*1*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < amphi_editing[j].position) {
                                                    if (modified_alignment_lines[amphi_order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[amphi_order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[amphi_order*2 + 1].insert(k, "*1*");
                                            }
                                            int order = 0;
                                            bool id_found = false;
                                            while (not id_found and order < order_of_ids.size()) {
                                                if (order_of_ids[order].compare(human_editing[hjs[hjj]].trans_id) == 0) id_found = true;
                                                else order++;
                                            }
                                            if (human_editing[hjs[hjj]].aligned_position - 1 == aligned_fasta_lines[order*2 + 1].length()) modified_alignment_lines[order*2 + 1].insert(aligned_fasta_lines[order*2 + 1].length(), "*1*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < human_editing[hjs[hjj]].position) {
                                                    if (modified_alignment_lines[order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[order*2 + 1].insert(k, "*1*");
                                            }
                                        }
                                        else if (amphi_editing[j].edi_aa == human_editing[hjs[hjj]].edi_aa) {
                                            class_2_out << amphi_editing[j].line << '\t' << amphi_editing[j].ori_aa << amphi_editing[j].position << amphi_editing[j].edi_aa << endl;
                                            class_2_out << human_editing[hjs[hjj]].line << endl;
                                            if (amphi_editing[j].aligned_position - 1 == aligned_fasta_lines[amphi_order*2 + 1].length()) modified_alignment_lines[amphi_order*2 + 1].insert(aligned_fasta_lines[amphi_order*2 + 1].length(), "*2*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < amphi_editing[j].position) {
                                                    if (modified_alignment_lines[amphi_order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[amphi_order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[amphi_order*2 + 1].insert(k, "*2*");
                                            }
                                            int order = 0;
                                            bool id_found = false;
                                            while (not id_found and order < order_of_ids.size()) {
                                                if (order_of_ids[order].compare(human_editing[hjs[hjj]].trans_id) == 0) id_found = true;
                                                else order++;
                                            }
                                            if (human_editing[hjs[hjj]].aligned_position - 1 == aligned_fasta_lines[order*2 + 1].length()) modified_alignment_lines[order*2 + 1].insert(aligned_fasta_lines[order*2 + 1].length(), "*2*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < human_editing[hjs[hjj]].position) {
                                                    if (modified_alignment_lines[order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[order*2 + 1].insert(k, "*2*");
                                            }
                                        }
                                        else if (amphi_editing[j].ori_aa == human_editing[hjs[hjj]].ori_aa) {
                                            class_3_out << amphi_editing[j].line << '\t' << amphi_editing[j].ori_aa << amphi_editing[j].position << amphi_editing[j].edi_aa << endl;
                                            class_3_out << human_editing[hjs[hjj]].line << endl;
                                            if (amphi_editing[j].aligned_position - 1 == aligned_fasta_lines[amphi_order*2 + 1].length()) modified_alignment_lines[amphi_order*2 + 1].insert(aligned_fasta_lines[amphi_order*2 + 1].length(), "*3*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < amphi_editing[j].position) {
                                                    if (modified_alignment_lines[amphi_order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[amphi_order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[amphi_order*2 + 1].insert(k, "*3*");
                                            }
                                            int order = 0;
                                            bool id_found = false;
                                            while (not id_found and order < order_of_ids.size()) {
                                                if (order_of_ids[order].compare(human_editing[hjs[hjj]].trans_id) == 0) id_found = true;
                                                else order++;
                                            }
                                            if (human_editing[hjs[hjj]].aligned_position - 1 == aligned_fasta_lines[order*2 + 1].length()) modified_alignment_lines[order*2 + 1].insert(aligned_fasta_lines[order*2 + 1].length(), "*3*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < human_editing[hjs[hjj]].position) {
                                                    if (modified_alignment_lines[order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[order*2 + 1].insert(k, "*3*");
                                            }
                                        }
                                        else {
                                            class_4_out << amphi_editing[j].line << '\t' << amphi_editing[j].ori_aa << amphi_editing[j].position << amphi_editing[j].edi_aa << endl;
                                            class_4_out << human_editing[hjs[hjj]].line << endl;
                                            if (amphi_editing[j].aligned_position - 1 == aligned_fasta_lines[amphi_order*2 + 1].length()) modified_alignment_lines[amphi_order*2 + 1].insert(aligned_fasta_lines[amphi_order*2 + 1].length(), "*4*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < amphi_editing[j].position) {
                                                    if (modified_alignment_lines[amphi_order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[amphi_order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[amphi_order*2 + 1].insert(k, "*4*");
                                            }
                                            int order = 0;
                                            bool id_found = false;
                                            while (not id_found and order < order_of_ids.size()) {
                                                if (order_of_ids[order].compare(human_editing[hjs[hjj]].trans_id) == 0) id_found = true;
                                                else order++;
                                            }
                                            if (human_editing[hjs[hjj]].aligned_position - 1 == aligned_fasta_lines[order*2 + 1].length()) modified_alignment_lines[order*2 + 1].insert(aligned_fasta_lines[order*2 + 1].length(), "*4*");
                                            else {
                                                int k = 0;
                                                int count = 0;
                                                while (count < human_editing[hjs[hjj]].position) {
                                                    if (modified_alignment_lines[order*2 + 1][k] == '*') k += 3;
                                                    else if (modified_alignment_lines[order*2 + 1][k++] != '-') count++;
                                                }
                                                modified_alignment_lines[order*2 + 1].insert(k, "*4*");
                                            }
                                        }
                                    }
                                }
                            }
                            for (int k = 1; k < modified_alignment_lines.size(); k += 2) {
                                if (modified_alignment_lines[k].length() != aligned_fasta_lines[k].length()) {
                                    aligned_out << modified_alignment_lines[k - 1] << endl;
                                    aligned_out << modified_alignment_lines[k] << endl;
                                }
                            }
                        }
                    }
                }
            }
            /*aligned_out << endl;
            class_0_out << endl;
            class_1_out << endl;
            class_2_out << endl;
            class_3_out << endl;
            class_4_out << endl;*/
            human_editing.clear();
            amphi_editing.clear();
            amphi_trans_ids.clear();
            human_genes.clear();
        }
        else if (var_line[0] == 'c') {
            vector<string> human_parsed = split(var_line, '\t');
            if (human_parsed[9].compare("exonic") == 0 and human_parsed[12].compare("UNKNOWN") != 0) {
                editing_info n_editing;
                n_editing.line = var_line;
                n_editing.gene_id = human_parsed[10];
                vector<string> editing_parsed = split(human_parsed[12], ':');
                n_editing.trans_id = editing_parsed[1];
                bool gene_found = false;
                int i = 0;
                while (not gene_found and i < human_genes.size()) {
                    gene_found = human_genes[i++].compare(n_editing.trans_id) == 0;
                }
                if (not gene_found) human_genes.push_back(n_editing.trans_id);
                vector<string> aa_parsed = split(editing_parsed[4], ',');
                string aa_info = aa_parsed[0];
                n_editing.edi_aa = aa_info[aa_info.length() - 1];
                aa_info.erase(aa_info.length() - 1, 1);
                n_editing.ori_aa = aa_info[2];
                aa_info.erase(0, 3);
                n_editing.position = atoi(aa_info.c_str());
                human_editing.push_back(n_editing);
            }
        }
        else {
            vector<string> amphi_parsed = split(var_line, '\t');
            vector<string> amphi_trans = split(amphi_parsed[12], ':');
            for (int t = 0; t < amphi_trans.size(); t++) {
                editing_info n_editing;
                n_editing.line = var_line;
                n_editing.gene_id = amphi_parsed[11];
                n_editing.trans_id = amphi_trans[t];
                int scaf_pos = atoi(amphi_parsed[1].c_str());
                map<string, vector<gtf_info> >::iterator gtf_it = amphi_gtf_info.find(n_editing.trans_id);
                if (gtf_it == amphi_gtf_info.end()) {
                    cerr << "Something went wrong, Unable to find gtf info on transcript " << n_editing.trans_id << endl; 
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
                        map<string, string>::iterator fas_it = amphi_trans_seqs.find(n_editing.trans_id);
                        if (fas_it == amphi_trans_seqs.end()) {
                            cerr << "Something went wrong, Unable to find fasta info on transcript " << n_editing.trans_id << endl; 
                        }
                        else {
                            bool wrong_nucleotide = false;
                            n_editing.position = trans_pos/3 + 1;
                            char first, second, third;
                            if (trans_pos%3 == 0) {
                                first = fas_it->second[trans_pos];
                                second = fas_it->second[trans_pos + 1];
                                third = fas_it->second[trans_pos + 2];
                                if (first != 'A' and first != 'a') {
                                    cerr << "Warning: Something went wrong, the nucleaotide in  " << n_editing.trans_id << " at position " << trans_pos << " is a " << first << ", not an A!" << endl; 
                                    wrong_nucleotide = true;
                                }
                                n_editing.ori_aa = translate_codon(first, second, third, genetic_code_table);
                                n_editing.edi_aa = translate_codon('G', second, third, genetic_code_table);
                            }
                            else if (trans_pos%3 == 1) {
                                first = fas_it->second[trans_pos - 1];
                                second = fas_it->second[trans_pos];
                                third = fas_it->second[trans_pos + 1];
                                if (second != 'A' and second != 'a') {
                                    cerr << "Warning: Something went wrong, the nucleaotide in  " << n_editing.trans_id << " at position " << trans_pos << " is a " << second << ", not an A!" << endl; 
                                    wrong_nucleotide = true;
                                }
                                n_editing.ori_aa = translate_codon(first, second, third, genetic_code_table);
                                n_editing.edi_aa = translate_codon(first, 'G', third, genetic_code_table);
                            }
                            else {
                                first = fas_it->second[trans_pos - 2];
                                second = fas_it->second[trans_pos - 1];
                                third = fas_it->second[trans_pos];
                                if (third != 'A' and third != 'a') {
                                    cerr << "Warning: Something went wrong, the nucleaotide in  " << n_editing.trans_id << " at position " << trans_pos << " is a " << third << ", not an A!" << endl; 
                                    wrong_nucleotide = true;
                                }
                                n_editing.ori_aa = translate_codon(first, second, third, genetic_code_table);
                                n_editing.edi_aa = translate_codon(first, second, 'G', genetic_code_table);
                            }
                            if (not wrong_nucleotide) {
                                amphi_editing.push_back(n_editing);
                                bool trans_found = false;
                                int i = 0;
                                while (not trans_found and i < amphi_trans_ids.size()) {
                                    trans_found = amphi_trans_ids[i++].compare(n_editing.trans_id) == 0;
                                }
                                if (not trans_found) amphi_trans_ids.push_back(n_editing.trans_id);
                            }
                        }
                    }
                }
            }
        }
    }
    var_in.close();
    class_0_out.close();
    class_1_out.close();
    class_2_out.close();
    class_3_out.close();
    class_4_out.close();
 //   if (aaRange != 0) class_5_out.close();
    aligned_out.close();
}