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
            cout << "Usage: " << argv[0] << "filtered_variants_amphioxus RNAediting_homo < homology_map" << endl;
            return 0;
        }
    }
    if (argc != 3) {
            cout << "Usage: " << argv[0] << "filtered_variants_amphioxus RNAediting_homo < homology_map" << endl;
            return 1;
    }
    const string var_in_path(argv[1]);
    const string human_in_path(argv[2]);
    bool debug = false;
    ifstream var_in;
    ifstream human_in;
    var_in.open(var_in_path.c_str());
    if (not var_in.is_open()) {
        cerr << "could not open filtered variants file" << endl;
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
    string amphi_gene;
    vector<string> homologs_found;
    map<string, vector<vector<int> >, comparer> homolog_editing;
    while (cin >> amphi_gene) {
        map<string, vector< vector<int> > >::iterator it = homolog_editing.insert(pair<string, vector< vector<int> > > (amphi_gene, vector<vector<int> > ())).first;
        int size1;
        cin >> size1;
        for (int i = 0; i < size1; i++) {
            it->second.push_back(vector<int> ());
            int size2;
            cin >> size2;
            for (int j = 0; j < size2; j++) {
                int line_n;
                cin >> line_n;
                it->second[i].push_back(line_n);
            }
        }
    }
    string prev_gene_id = "";
    string var_line;
    bool eof = false;
    if (not getline(var_in, var_line)) {
        eof = true;
    }
    while (not eof) {
        vector<string> parsed = split(var_line, '\t');
        string gene_id = parsed[11];
        if (gene_id.compare("unknown_gene") != 0 and gene_id.compare(prev_gene_id) != 0) {
            prev_gene_id = gene_id;
            map<string, vector< vector<int> > >::iterator it;
            it = homolog_editing.find(gene_id);
            if (it != homolog_editing.end()) {
                while (gene_id == prev_gene_id and not eof) {
                    cout << var_line << endl;
                    if (not getline(var_in, var_line)) eof = true;
                    else {
                        parsed = split(var_line, '\t');
                        gene_id = parsed[11];
                    }
                }
                for (int i = 0; i < it->second.size(); i++) {
                    for (int j = 0; j < it->second[i].size(); j++) {
                        cout << human_editing[it->second[i][j]] << endl;
                    }
                }
                cout << endl;
            }
            else {
                if (not getline(var_in, var_line)) eof = true;
            }
        }
        else {
            if (not getline(var_in, var_line)) eof = true;
        }
    }
}    
/*    string prev_gene_id = "";
    string var_line;
    bool eof = false;
    if (not getline(var_in, var_line)) {
        eof = true;
    }
    bool eof2 = false;
    string blast_line;
    if (not getline(blast_in, blast_line)) eof2 = true;
    while (not eof) {
        if (var_line.compare("Sc0000000	3320948	.	T	C	225.009	.	DP=310;VDB=0.954996;SGB=-0.693147;RPB=0.865703;MQB=0.000309409;MQSB=0.790887;BQB=0.998136;MQ0F=0;AF1=0.5;AC1=1;DP4=110,40,75,22;MQ=38;FQ=225.007;PV4=0.548757,0.331929,0.00107179,1	GT:PL	0/1:255,0,255	-	BL02686	BL02686_evm7") == 0) debug = true;
        if (debug) cout << "debug: " << endl;
        if (debug) cout << "debug: " << "analyzing line: " << var_line << endl;
        vector<string> parsed = split(var_line, '\t');
        string gene_id = parsed[11];
        if (debug) cout << "debug: " << gene_id << endl;
        if (gene_id.compare("unknown_gene") != 0 and gene_id.compare(prev_gene_id) != 0) {
            prev_gene_id = gene_id;
            bool homolog_found = false;
            vector<string> human_genes_checked;
            while (not homolog_found and not eof2) {
                if (debug) cout << "debug: " << "analyzing homology for: " << blast_line << endl;
                if (blast_line[0] != '#') {
                    vector<string> blast_parsed = split(blast_line, '\t');
                    vector<string> trans_id_parsed = split(blast_parsed[0], '_');
                    
                    if (debug) cout << "debug: " << trans_id_parsed[0];
                    
                    if (trans_id_parsed[0].compare(gene_id) == 0) {
                        if (debug) cout << "debug: " << "homolog found!" << endl;
                        if (debug) cout << "debug: " << blast_line << endl;
                        vector<string> human_id_parsed = split(blast_parsed[1], '|');
                        string human_gene = human_id_parsed[2];
                        bool humid_checked = false;
                        int i = 0;
                        while (not humid_checked and i < human_genes_checked.size()) {
                            if (human_genes_checked[i].compare(human_gene) == 0) humid_checked = true;
                            i++;
                        }
                        if (not humid_checked) {
                            human_genes_checked.push_back(human_gene);
                            human_in.open(human_in_path.c_str());
                            if (not human_in.is_open()) {
                                cerr << "could not open fasta file" << endl;
                                return 1;
                            }
                            string human_line;
                            while (getline(human_in, human_line)) {
                                vector<string> human_parsed = split (human_line,'\t');
                                if (human_parsed[9].compare("exonic") == 0 and human_parsed[10].compare(human_gene) == 0
                                or human_parsed[14].compare("exonic") == 0 and human_parsed[15].compare(human_gene) == 0
                                or human_parsed[18].compare("exonic") == 0 and human_parsed[19].compare(human_gene) == 0) {
                                    if (not homolog_found) {
                                        homolog_found = true;
                                        while (gene_id == prev_gene_id and not eof) {
                                            cout << var_line << endl;
                                            if (not getline(var_in, var_line)) eof = true;
                                            else {
                                                parsed = split(var_line, '\t');
                                                gene_id = parsed[11];
                                            }
                                        }
                                    }
                                    cout << human_line << endl;
                                }
                            }
                            human_in.close();
                        }
                        if (not getline(blast_in, blast_line)) eof2 = true;
                        if (debug) cout << "debug: " << "end of homology analysis" << endl;
                    }
                    else {
                        if (not getline(blast_in, blast_line)) eof2 = true;
                    }
                        
                }
                else {
                    if (not getline(blast_in, blast_line)) eof2 = true;
                }
                    
            }
            if (not homolog_found) {
                if (debug) cout << "debug: " << "searching again" << endl;
                blast_in.close();
                blast_in.open(blast_in_path.c_str());
                if (not blast_in.is_open()) {
                    cerr << "could not open blast results file" << endl;
                    return 1;
                }
                eof2 = false;
                if (not getline(blast_in, blast_line)) eof2 = true;
                while (not homolog_found and not eof2) {
                    if (blast_line[0] != '#') {
                        vector<string> blast_parsed = split(blast_line, '\t');
                        vector<string> trans_id_parsed = split(blast_parsed[0], '_');
                        //if (debug) cout << "debug: " << trans_id_parsed[0] << endl;
                        if (trans_id_parsed[0].compare(gene_id) == 0) {
                            if (debug) cout << "debug: " << "homolog found! 2" << endl;
                            if (debug) cout << "debug: " << blast_line << endl;

                            vector<string> human_id_parsed = split(blast_parsed[1], '|');
                            string human_gene = human_id_parsed[2];
                            bool humid_checked = false;
                            int i = 0;
                            while (not humid_checked and i < human_genes_checked.size()) {
                                if (human_genes_checked[i].compare(human_gene) == 0) humid_checked = true;
                                i++;
                            }
                            if (not humid_checked) {
                                human_genes_checked.push_back(human_gene);
                                human_in.open(human_in_path.c_str());
                                if (not human_in.is_open()) {
                                    cerr << "could not open fasta file" << endl;
                                    return 1;
                                }
                                string human_line;
                                while (getline(human_in, human_line)) {
                                    vector<string> human_parsed = split (human_line,'\t');
                                    if (human_parsed[9].compare("exonic") == 0 and human_parsed[10].compare(human_gene) == 0
                                    or human_parsed[14].compare("exonic") == 0 and human_parsed[15].compare(human_gene) == 0
                                    or human_parsed[18].compare("exonic") == 0 and human_parsed[19].compare(human_gene) == 0) {
                                        if (not homolog_found) {
                                            homolog_found = true;
                                            while (gene_id == prev_gene_id and not eof) {
                                                cout << var_line << endl;
                                                if (not getline(var_in, var_line)) eof = true;
                                                else {
                                                    parsed = split(var_line, '\t');
                                                    gene_id = parsed[11];
                                                }
                                            }
                                        }
                                        cout << human_line << endl;
                                    }
                                }
                                human_in.close();
                            }
                            if (not getline(blast_in, blast_line)) eof2 = true;
                        }
                        else {
                            if (not getline(blast_in, blast_line)) eof2 = true;
                        }
                            
                    }
                    else {
                        if (not getline(blast_in, blast_line)) eof2 = true;
                    }
                }
            }
            if (homolog_found) {
                if (debug) cout << "debug: " << "case homolog found" << endl;
                if (not eof2) {
                    if (debug) cout << "debug: " << blast_line << endl;
                    if (debug) cout << "debug: " << "NOT EOF2" << endl;
                    vector<string> blast_parsed = split(blast_line, '\t');
                    if (debug) cout << "debug: " << "blast_parsed" << endl;
                    vector<string> trans_id_parsed = split(blast_parsed[0], '_');
                    if (debug) cout << "debug: " << "trans_id_parsed" << endl;
                    if (debug) cout << "debug: " << blast_parsed[0] << endl;
                    if (debug) cout << "debug: " << blast_parsed[1] << endl;
                    vector<string> human_id_parsed = split(blast_parsed[1], '|');
                    if (debug) cout << "debug: " << "human_id_parsed" << endl;
                    string human_gene = human_id_parsed[2];
                    if (debug) cout << "debug: " << "human_gene" << endl;
                    int n_homologs = 1;
                    while (trans_id_parsed[0].compare(gene_id) == 0 and not eof2 and n_homologs < 5) {
                        if (debug) cout << "debug: " << blast_line << endl;
                        bool humid_checked = false;
                        int i = 0;
                        while (not humid_checked and i < human_genes_checked.size()) {
                            if (human_genes_checked[i].compare(human_gene) == 0) humid_checked = true;
                            i++;
                        }
                        if (not humid_checked) {
                            human_genes_checked.push_back(human_gene);
                            human_in.open(human_in_path.c_str());
                            if (not human_in.is_open()) {
                                cerr << "could not open fasta file" << endl;
                                return 1;
                            }
                            string human_line;
                            int prev_n_homologs = n_homologs;
                            while (getline(human_in, human_line)) {
                                if (debug) cout << "debug: " << human_line << endl;
                                vector<string> human_parsed = split (human_line,'\t');
                                if (human_parsed[9].compare("exonic") == 0 and human_parsed[10].compare(human_gene) == 0
                                or human_parsed[14].compare("exonic") == 0 and human_parsed[15].compare(human_gene) == 0
                                or human_parsed[18].compare("exonic") == 0 and human_parsed[19].compare(human_gene) == 0) {
                                    cout << human_line << endl;
                                    if (n_homologs == prev_n_homologs) {
                                        n_homologs++;
                                    }
                                }
                            }
                            if (debug) cout << "debug: " << "end human loop" << endl;
                            human_in.close();
                        }
                        if (not getline(blast_in, blast_line)) eof2 = true;
                        else {
                            blast_parsed = split(blast_line, '\t');
                            trans_id_parsed = split(blast_parsed[0], '_');
                            human_id_parsed = split(blast_parsed[1], '|');
                            human_gene = human_id_parsed[2];
                        }
                    }
                    if (debug) cout << "debug: " << "end blastloop" << endl;                    
                }
                cout << endl;
            }
            else {
                if (not getline(var_in, var_line)) eof = true;
            }
        }
        else {
            if (not getline(var_in, var_line)) {
                eof = true;
            }
        }
    }
    var_in.close();
    blast_in.close();
}*/