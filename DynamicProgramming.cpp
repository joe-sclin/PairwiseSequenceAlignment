#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>

using namespace std;

signed int Match = 1;       // Initialise scoring variable using standard value
signed int Mismatch = -1;   // Match and mismatch only applied to DNA / RNA alignment
signed int Gap = -2;
vector<vector<int>> aa_matrix;  // Only applied to protein alignment
map<string , int> aa_mapping;
const int HOME = 0;        // Global const variable for direction in scoring matrix (0 = binary 00)
const int DIAGONAL = 1;    // (binary 01)
const int HORIZONTAL = 2;  // (binary 10)
const int VERTICAL = 4;    // (binary 100)


signed int base_score(char a, char b = '\0') {  // Score function for DNA / RNA
    return (a == '\0' && b == '\0') ? Gap : (a == b) ? Match : Mismatch;
}
signed int aa_score(char a, char b) {       // Score function for protein (Refer to mapping on matrix --> Slower than static value in DNA / RNA)
    if (a == '\0' && b == '\0') {
        return aa_matrix[aa_mapping.at("*")][aa_mapping.at("*")];
    } else {
        string str_a(1, a);
        string str_b(1, b);
        return aa_matrix[aa_mapping.at(str_a)][aa_mapping.at(str_b)];
    }
}

void BLOSUM_matrix() {   // Modify BLOSUM_matrix for Match, Mismatch and Gap score
    ifstream BLOSUM_file("BLOSUM62.csv");
    vector<vector<int>> BLOSUM62;       // Vector can only store matrix value
    map<string , int> aa_map;
    string currentline;
    while (getline(BLOSUM_file, currentline)) {
        stringstream ss(currentline);
        // Skip the first row
        if (currentline[0] == ',') {    // First row (List of amino acid)
            int cols = -1; // First row and column are skipped (space is mapped to -1)
            string aminoacid;
            while (getline(ss, aminoacid, ',')) {
                aa_map[aminoacid] = cols;
                cols++;
            }
            continue;
        }
        vector<int> currentrow;
        string value;
        getline(ss, value, ',');
        while (getline(ss, value, ',')) {
            currentrow.push_back(stoi(value));
        }
        BLOSUM62.push_back(currentrow);
    }
    BLOSUM_file.close();
    aa_matrix = BLOSUM62;
    aa_mapping = aa_map;
}

// Traceback function to generate aligned sequences by dynamic programming approach (Not storing the traceback path, only storing aligned part of sequences)
pair<string, string> traceback(unsigned int i, unsigned int j, vector<vector<int>>& T, string& x, string& y, char GAP = '-') {
    string Ax, Ay;  // Initialise Aligned seq x (Ax) and seq y (Ay)
    while (T[i][j] != HOME) {
        int trace = T[i][j];
        if (trace & DIAGONAL) {     // Appending base for both sequences
            i--; Ax += x[i];
            j--; Ay += y[j];
            continue;
        }
        if (trace & HORIZONTAL) {   // Only appending base for sequence y (x = gap)
            Ax += GAP;
            j--; Ay += y[j];
            continue;
        }
        if (trace & VERTICAL) {   // Only appending base for sequence x (y = gap)
            i--;
            Ax += x[i];
            Ay += GAP;
        }
    }
    // As traceback is started from the last base of sequences -> reverse order by rbegin and rend
    return make_pair(string(Ax.rbegin(), Ax.rend()), string(Ay.rbegin(), Ay.rend()));   // Pair of aligned sequences x and y
}

// Section: implementation of 4 alignment variants
pair<string, string> global(string& x, string& y, int(*score)(char, char)) {     // Needleman-Wunsch algorithm
    unsigned int m = x.length();
    unsigned int n = y.length();
    vector<vector<int>> Sc(m+1, vector<int>(n+1)); // scores matrix
    vector<vector<int>> Trace(m+1, vector<int>(n+1)); // traceback matrix
    Sc[0][0] = HOME;        // Fill starting point
    for (int j = 1; j <= n; j++) {      // Fill first row (Horizontal topmost)
        Sc[0][j] = j * Gap;
        Trace[0][j] = HORIZONTAL;
    }
    for (int i = 1; i <= m; i++) {      // Fill first column (Vertical leftmost)
        Sc[i][0] = i * Gap;
        Trace[i][0] = VERTICAL;
    }
    for (int i = 1; i <= m; i++) {      // Fill remaining cells in matrix Sc
        char xi = x[i-1];
        for (int j = 1; j <= n; j++) {
            char yj = y[j-1];
            int d = Sc[i-1][j-1] + (*score)(xi, yj);
            int h = Sc[i][j-1] + Gap;
            int v = Sc[i-1][j] + Gap;
            int best_dir = max({d, h, v});      // Select the direction with the best score (the highest value)
            Sc[i][j] = best_dir;
            Trace[i][j] = ((d==best_dir)*DIAGONAL) + ((h==best_dir)*HORIZONTAL) + ((v==best_dir)*VERTICAL);
        }
    }
    int finalscore = Sc[m][n];
    cout << "Alignment score: " << finalscore << endl;
    pair<string, string> alignedseq = traceback(m, n, Trace, x, y);
    return alignedseq;
}

pair<string, string> local(string& x, string& y, int(*score)(char, char)) {     // Smith-Waterman algorithm
    unsigned int m = x.length();
    unsigned int n = y.length();
    vector<vector<int>> Sc(m+1, vector<int>(n+1));
    vector<vector<int>> Trace(m+1, vector<int>(n+1));
    Sc[0][0] = HOME;        // Fill with 0
    for (int j = 1; j <= n; j++) {      // Fill first row with 0
        Trace[0][j] = HOME;
    }
    for (int i = 1; i <= m; i++) {      // Fill first column with 0
        Trace[i][0] = HOME;
    }
    // Initialise tuple (alignment score, aligned sequence x, aligned sequence y)
    tuple<int, int, int> best = make_tuple(-1, -1, -1);
    for (int i = 1; i <= m; i++) {      // Fill remaining cells in matrix Sc
        char xi = x[i-1];
        for (int j = 1; j <= n; j++) {
            char yj = y[j-1];
            int d = Sc[i-1][j-1] + (*score)(xi, yj);
            int h = Sc[i][j-1] + Gap;
            int v = Sc[i-1][j] + Gap;
            int best_dir = max({0, d, h, v});      // 0 is also included (Not in global alignment)
            Sc[i][j] = best_dir;
            Trace[i][j] = ((d==best_dir)*DIAGONAL) + ((h==best_dir)*HORIZONTAL) + ((v==best_dir)*VERTICAL);
            if (Sc[i][j] > get<0>(best)){   // update the best aligned region
                get<0>(best) = Sc[i][j];
                get<1>(best) = i;
                get<2>(best) = j;
            }
        }
    }
    int finalscore = get<0>(best);
    m = get<1>(best);       // Local alignment -> length of aligned part may not equal to full sequence
    n = get<2>(best);
    cout << "Alignment score: " << finalscore << endl;
    pair<string, string> alignedseq = traceback(m, n, Trace, x, y);
    return alignedseq;
}

pair<string, string> semiglobal(string& x, string& y, int(*score)(char, char)) {     // Modified global alignment (No gap penalty at the end of either sequences)
    unsigned int m = x.length();
    unsigned int n = y.length();
    vector<vector<int>> Sc(m+1, vector<int>(n+1));
    vector<vector<int>> Trace(m+1, vector<int>(n+1));
    Sc[0][0] = HOME;
    Trace[0][0] = HOME;
    for (int j = 1; j <= n; j++) {
        Sc[0][j] = j * Gap;
        Trace[0][j] = HOME;
    }
    for (int i = 1; i <= m; i++) {
        Sc[i][0] = i * Gap;
        Trace[i][0] = VERTICAL;
    }
    tuple<int, int, int> best = make_tuple(-1, -1, -1);
    for (int i = 1; i <= m; i++) {
        char xi = x[i-1];
        for (int j = 1; j <= n; j++) {
            char yj = y[j-1];
            int d = Sc[i-1][j-1] + (*score)(xi, yj);
            int h = Sc[i][j-1] + Gap;
            int v = Sc[i-1][j] + Gap;
            int best_dir = max({ d, h, v});
            Sc[i][j] = best_dir;
            Trace[i][j] = ((d==best_dir)*DIAGONAL) + ((h==best_dir)*HORIZONTAL) + ((v==best_dir)*VERTICAL);
            if (Sc[m][j] > get<0>(best)){
                get<0>(best) = Sc[m][j];
                get<1>(best) = i;
                get<2>(best) = j;
            }
        }
    }
    int finalscore = get<0>(best);
    cout << "Alignment score: " << finalscore << endl;
    pair<string, string> alignedseq = traceback(m, n, Trace, x, y);
    return alignedseq;
}

pair<string, string> overlap(string& x, string& y, int(*score)(char, char)) {     // Modified local alignment (No gap penalty at the end of either sequences)
    unsigned int m = x.length();
    unsigned int n = y.length();
    vector<vector<int>> Sc(m + 1, vector<int>(n + 1));
    vector<vector<int>> Trace(m + 1, vector<int>(n + 1));
    Trace[0][0] = HOME;
    for (int j = 1; j <= n; j++) {
        Trace[0][j] = HOME;
    }
    for (int i = 1; i <= m; i++) {
        Trace[i][0] = HOME;
    }
    tuple<int, int, int> best = make_tuple(-1, -1, -1);
    for (int i = 1; i <= m; i++) {
        char xi = x[i - 1];
        for (int j = 1; j <= n; j++) {
            char yj = y[j - 1];
            int d = Sc[i - 1][j - 1] + (*score)(xi, yj);
            int h = Sc[i][j - 1] + Gap;
            int v = Sc[i - 1][j] + Gap;
            if (i == 0 || j == 0){
                int best_dir = max({0, d, h, v});
                Sc[i][j] = best_dir;
                Trace[i][j] = ((d == best_dir) * DIAGONAL) + ((h == best_dir) * HORIZONTAL) + ((v == best_dir) * VERTICAL);
            } else {
                int best_dir = max({d, h, v});
                Sc[i][j] = best_dir;
                Trace[i][j] = ((d == best_dir) * DIAGONAL) + ((h == best_dir) * HORIZONTAL) + ((v == best_dir) * VERTICAL);
            }
            if (Sc[i][j] > get<0>(best)) {
                get<0>(best) = Sc[i][j];
                get<1>(best) = i;
                get<2>(best) = j;
            }
        }
    }
    int finalscore = get<0>(best);
    int i = get<1>(best);
    int j = get<2>(best);
    cout << "Alignment score: " << finalscore << endl;
    pair<string, string> alignedseq = traceback(i, j, Trace, x, y);
    return alignedseq;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {    // Check number of input argument
        cerr << "Error: Format of input: " << argv[0] << " sequence.fasta" << endl;
        exit(1);
    } else {
        // Section: fasta file checking
        ifstream fasta(argv[1]);

        if (!fasta.is_open()) {     // Check whether fasta files can open
            cerr << "Error: could not open fasta file" << endl;
            exit(1);
        }
        vector<string> sequences;   // Temp object to store sequences read
        string CurrentLine;

        while (getline(fasta, CurrentLine)) {
            if (CurrentLine.empty()) {
                continue; // Skip empty lines between sequences
            }
            if (CurrentLine[0] ==
                '>') {    // Name of sequences are omitted in alignment (But > symbol represented start of another sequence)
                sequences.emplace_back();   // Create an empty new element in vector for appending current line
            } else if (sequences.empty()) {
                cerr << "Error: No sequence present in fasta file / Incorrect format of fasta file" << endl;
            } else {
                sequences.back() += CurrentLine; // Append the current line to the current sequence
            }
        }
        if (sequences.size() < 2) {
            cerr << "Error: At least 2 sequences are required in input file" << endl;
            return 1;
        } else if (sequences.size() > 2) {
            cout << "Reminder: Only pairwise alignment between first two sequences will be performed." << endl;
        }
        string seq1 = sequences[0];
        string seq2 = sequences[1];
        // convert all characters in seq1 and seq 2 to uppercase (Case insensitive alignment program)
        transform(seq1.begin(), seq1.end(), seq1.begin(), ::toupper);
        transform(seq2.begin(), seq2.end(), seq2.begin(), ::toupper);

        // Section: input sequences checking
        string choice;
        cout << "Which type of alignment do you want (DNA / RNA / protein)" << endl;
        cin >> choice;
        if (choice != "DNA" && choice != "RNA" && choice != "protein") {
            cerr << "Error: Only DNA / RNA / protein is allowed." << endl;
            exit(1);
        } else if (choice == "DNA") {   // check that all characters are A/T/C/G
            for (char c: seq1) {
                if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
                    cerr << "Error: Invalid character in Sequence 1, only A / C / G / T are allowed." << endl;
                    exit(1);
                }
            }
            for (char c: seq2) {
                if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
                    cerr << "Error: Invalid character in Sequence 2, only A / C / G / T are allowed." << endl;
                    exit(1);
                }
            }
        } else if (choice == "RNA") {   // check that all characters are A/U/C/G
            for (char c: seq1) {
                if (c != 'A' && c != 'U' && c != 'C' && c != 'G') {
                    cerr << "Error: Invalid character in Sequence 1, only A / C / G / U are allowed." << endl;
                    exit(1);
                }
            }
            for (char c: seq2) {
                if (c != 'A' && c != 'U' && c != 'C' && c != 'G') {
                    cerr << "Error: Invalid character in Sequence 2, only A / C / G / U are allowed." << endl;
                    exit(1);
                }
            }
        } else {    // protein
            const string aa_list = "ARNDCQEGHIKLMFPSTWYVXBZX";
            for (char c: seq1) {
                if (aa_list.find(c) == string::npos) {
                    cerr << "Error: Invalid character in Sequence 1, only general amino acid codon are allowed." << endl;
                }
            }
            for (char c: seq2) {
                if (aa_list.find(c) == string::npos) {
                    cerr << "Error: Invalid character in Sequence 1, only general amino acid codon are allowed." << endl;
                }
            }
        }

        // Section: Select different score functions for DNA/RNA or protein
        int (*selected_score)(char, char);  // Using function pointer to select different score function for DNA/RNA or protein alignment
        // Section: change scores setting
        if (choice == "DNA" || choice == "RNA") {
            selected_score = base_score;
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');    // Clear '\n' from previous question
            char change;
            cout << "Do you want to change default score setting (Y or N)? (Match: " << Match << ", Mismatch: "
                 << Mismatch
                 << ", Gap: " << Gap << endl;
            cin >> change;
            if (change != 'Y' && change != 'N') {
                cerr << "Error: Only Y or N is allowed." << endl;
                exit(1);
            } else if (change == 'Y') {
                cout << "Please enter the new Match, Mismatch, Gap score:" << endl;
                cin >> Match >> Mismatch >> Gap;
                if (cin.fail()) {
                    cerr << "Error: Only integer scores are allowed." << endl;
                    exit(1);
                }
            }
        } else {
            BLOSUM_matrix();
            selected_score = aa_score;
        }
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');    // Clear '\n' from previous question
        // Section: selection of alignment algorithm
        string alignment;
        cout << "Which alignment would you want to compute for input sequences: (global, local, semi-global or overlap)" << endl;
        getline(cin, alignment);
        pair<string, string> res;
        if (alignment == "global") {
            res = global(seq1, seq2, *selected_score);  // input function pointer of selected score function to alignment algorithm
        } else if (alignment == "local") {
            res = local(seq1, seq2, *selected_score);
        } else if (alignment == "semi-global") {
            res = semiglobal(seq1, seq2, *selected_score);
        } else if (alignment == "overlap") {
            res = overlap(seq1, seq2, *selected_score);
        } else {
            cerr << "Error: Invalid input alignment method, only global, local, semi-global or overlap are allowed."
                 << endl;
            exit(1);
        }
        cout << res.first << endl;
        cout << res.second << endl;
        std::string output_filename = string(argv[1]) + "_alignment_result_" + alignment + ".txt";
        std::ofstream output(output_filename);
        unsigned int SequenceLength = max(res.first.length(), res.second.length());
        int startpos = 0;       // Divide aligned region into 100bp per line
        int chunksize = 100;
        while (startpos < SequenceLength) {
            string CurrentChunk1 = res.first.substr(startpos, chunksize);
            string CurrentChunk2 = res.second.substr(startpos, chunksize);
            string MatchingLine(chunksize, ' ');
            for (int i = 0; i < min(CurrentChunk1.length(), CurrentChunk2.length()); i++) {
                if (CurrentChunk1[i] == CurrentChunk2[i]) {
                    MatchingLine[i] = '|';      // Only add matching line "|" when corresponding character is matched (Never have both gap in pairwise alignment)
                }
            }
            output << CurrentChunk1 << endl;
            output << MatchingLine << endl;
            output << CurrentChunk2 << endl;
            output << endl;
            startpos += chunksize;
        }
        output.close();
        return 0;
    }
}
