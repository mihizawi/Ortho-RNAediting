/*
 * This program creates a map using the results of a blast to relate the amphioxus genes with the lines of an RNA editing database that correspond to editing events in
 * the five best human homologs. The blast input file should be the amphiozus genes against the human genes and should have the standard blast output but with the lines of the 
 * same amphioxus gene grouped together and ordered by the blast score, having the best results for the amphioxus gene appear first (this can be achieved by sorting the lines by
 * the score and then sorting again with a stable algorithm by the amphixous gene). The human editing input file should have the format used in the REDIportal database, and
 * specifically for this program must have the gene id in the 11th column and informationn on whether the position is exonic on the 10th column (columns separated by tabs).
 * The output is a direct output of the map created by this program and is intended to be used as the input for the program cross_human_v2, together with the human editing file.
 * Note that the human editing file used here and in the cross_human_v2 programs must be the same. If another human input file is to be used, this program must be rerun to create
 * the corresponding map file. However, if the homology information (blast file) and the human editing file remain the same, then the map can be reused for multiple executions of
 * the cross_human_v2 program.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <stdlib.h>
using namespace std;

//split a string into a vector using a character as a delimitator, not original code
template<typename Out>
void split(const string &s, char delim, Out result) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

//split a string into a vector using a character as a delimitator, not original code
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

//a comparer for strings in function version as required for the map container
struct comparer
{
    public:
    bool operator()(const std::string x, const std::string y)
    {
         return x.compare(y) < 0;
    }
};

//the main function that will get called when the program is run
int main (int argc, char* argv[]) {
    //implementing -h option
    if (argc == 2) {
        string par(argv[1]);
        if (par == "-h") {
            cout << "Usage: " << argv[0] << "RNAediting_homo blast_amphi_homo" << endl;
            return 0;
        }
    }
    // checking if the minimum number of arguments is correct
    if (argc != 3) {
            cout << "Usage: " << argv[0] << "RNAediting_homo blast_amphi_homo" << endl;
            return 1;
    }
    // reading the paths of the input files from the two first arguments
    const string human_in_path(argv[1]);
    const string blast_in_path(argv[2]);
    // opening the human editing input file
    ifstream human_in;
    human_in.open(human_in_path.c_str());
    if (not human_in.is_open()) {
        cerr << "could not open human editing file" << endl;
        return 1;
    }
    //reading the human editing file and storing every line in a vector
    vector<string> human_editing;
    string human_line;
    while (getline(human_in, human_line)) {
        human_editing.push_back(human_line);
    }
    human_in.close();
    // declaring variables to read the blast lines
    string blast_line;
    string current_gene = "no_gene";
    // this vector will store for the current amphioxus gene, all the human homologs that contain an editing event (maximum 5)
    vector<string> homologs_found;
    // this map will store for each amphioxus gene (used as the key), a structure that stores all the numbers of the human editing file that represent
    // to editing events in the 5 best human homologs, grouped by each homologous gene
    map<string, vector<vector<int> >, comparer> homolog_editing;
    // opening the blast input file
    ifstream blast_in;
    blast_in.open(blast_in_path.c_str());
    if (not blast_in.is_open()) {
        cerr << "could not open blast results file" << endl;
        return 1;
    }
    //reading each blast line
    while (getline(blast_in, blast_line)) {
        //ignoring header or comment lines
        if (blast_line[0] != '#') {
            //parsing the blast line
            vector<string> blast_parsed = split(blast_line, '\t');
            //getting the amphioxus gene id from the transcript id
            vector<string> trans_id_parsed = split(blast_parsed[0], '_');
            string amphi_gene = trans_id_parsed[0];
            // checking if the amphioxus gene is different than that of the previous line
            if (current_gene.compare(amphi_gene) != 0) {
                // if the gene we read is different than the previous one, we have finished processing the previous one,  we retrieve the information
                // from the map and print it to the output
                map<string, vector< vector<int> > >::iterator it;
                it = homolog_editing.find(current_gene);
                if (it != homolog_editing.end()) {
                    int size1 = it->second.size();
                    // we first write the amphioxus gene and the number of different homolog genes with editing we have stored
                    cout << current_gene << " " << size1;
                    //for each homologous gene we write the number of editing lines stored for that gene and the line numbers
                    for (int i = 0; i < size1; i++) {
                        int size2 = it->second[i].size();
                        cout << " " << size2;
                        for (int j = 0; j < size2; j++) {
                            cout << " " << it->second[i][j];
                        }
                    }
                    cout << endl;
                }
                //after writing the stored information about the previous amphioxus gene, we update the current gene and we reset the list of human homologs
                homologs_found.clear();
                current_gene = amphi_gene;
            }
            //we check if we already have five homologs stored for the current amphioxus gene
            if (homologs_found.size() < 5) {
                //we get the human gene from the blast line
                vector<string> human_id_parsed = split(blast_parsed[1], '|');
                string human_gene = human_id_parsed[2];
                //we check if we already have stored information about this human gene
                bool seen_homolog = false;
                for (int i = 0; i < homologs_found.size(); i++) {
                    if (homologs_found[i].compare(human_gene) == 0) seen_homolog = true;
                }
                //if the human gene was not stored in association to the current amphioxus gene, we store the new information
                if (not seen_homolog) {
                    //we add the human gene to the vector storing the human homologs we have found
                    homologs_found.push_back(human_gene);
                    //we check if the map already has an entry for the current amphoxus gene
                    map<string, vector< vector<int> > >::iterator it;
                    it = homolog_editing.find(amphi_gene);
                    if (it == homolog_editing.end()) {
                        //if no entry is found for the current amphioxus gene, we create the whole data structure to be inserted for the amphioxus gene
                        vector<vector<int> > v;
                        v.push_back(vector<int> ());
                        //we search every human editing line from the input file that corresponds to the human gene we have found in the blast line
                        for (int i = 0; i < human_editing.size(); i++) {
                            vector<string> human_parsed = split(human_editing[i],'\t');
                            if (human_parsed[9].compare("exonic") == 0 and human_parsed[10].compare(human_gene) == 0) {
                                //we store the line number of the editing event
                                v[0].push_back(i);
                            }
                        }
                        //if we have found at least one line corresponding to the current human gene, we insert the structure for the amphioxus gene in the map
                        if (v[0].size() > 0) {
                            homolog_editing.insert(pair<string, vector< vector<int> > > (amphi_gene, v));
                        }
                    }
                    else {
                        //if there's already information about the amphioxus gene in the map, we only create the vector for the current human gene to be added to the structure
                        vector<int> v;
                        //we search every human editing line from the input file that corresponds to the human gene we have found in the blast line
                        for (int i = 0; i < human_editing.size(); i++) {
                            vector<string> human_parsed = split (human_editing[i],'\t');
                            if (human_parsed[9].compare("exonic") == 0 and human_parsed[10].compare(human_gene) == 0) {
                                v.push_back(i);
                            }
                        }
                        //if we have found at least one line corresponding to the current human gene, we add the vector to the structure stored in the map
                        if (v.size() > 0) {
                            it->second.push_back(v);
                        }
                    }
                }
            }
        }
    }
    blast_in.close();
}