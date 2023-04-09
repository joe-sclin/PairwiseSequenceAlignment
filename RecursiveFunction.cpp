#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

int Match = 1;
int Mismatch = -1;
int Gap = -2;
vector<vector<int>> Sc;

int AlignmentScore(string Seq1, string Seq2, int m, int n)      //function to find the best alignment score
{
    for (int i = 0; i <= m; i++)    //fill first row
        Sc[0][i] = Gap * i;

    for (int i = 0; i <= n; i++)    //fill first column
        Sc[i][0] = Gap * i;

    for (int i = 1; i <= m; i++) {      // filling remaining cells
        for (int j = 1; j <= n; j++){
            int score = (Seq1[i-1] == Seq2[j-1]) ? Match : Mismatch;
            Sc[i][j] = max({Sc[i-1][j-1] + score,
                               Sc[i-1][j] + Gap,
                               Sc[i][j-1] + Gap});
        }
    }
    return Sc[m][n];
}

//recursive function to find all possible branches
void OptimalAlignment(string Seq1, string Seq2, int m, int n, string curralignedX, string curralignedY)
{
    int i = m;
    int j = n;
    string alignmentX = curralignedX;
    string alignmentY = curralignedY;

    if(i == 0 && j == 0) {      // Complete alignment, output aligned result
        ofstream output("Recursive_all_optimal_alignment.txt", ios::app);       // Appending each alignment at the end
        output << alignmentX << endl;
        output << alignmentY << endl;
        output << std::endl;
    }
    else{   // Recursive part
            if (i == 0){ // Add free gap to AlignmentX
             OptimalAlignment(Seq1, Seq2, i, j - 1, "-" + alignmentX, Seq2[j-1] + alignmentY);
        } else if (j == 0){ // Add free gap to AlignmentY
            OptimalAlignment(Seq1, Seq2, i - 1, j, Seq1[i-1] + alignmentX, "-" + alignmentY);
        } else {    // Diagonal direction traceback
            int score1 = (Seq1[i-1] == Seq2[j-1]) ? Match : Mismatch;
            if (Sc[i][j] == Sc[i-1][j-1] + score1) {
                OptimalAlignment(Seq1, Seq2, i - 1, j - 1, Seq1[i-1] + alignmentX, Seq2[j-1] + alignmentY);
            }
            if (Sc[i][j] == Sc[i-1][j] + Gap) {     // Horizontal direction traceback (Add free gap to alignmentY)
                OptimalAlignment(Seq1, Seq2, i - 1, j, Seq1[i-1] + alignmentX, "-" + alignmentY);
            }
            if (Sc[i][j] == Sc[i][j-1] + Gap) {     // Vertical direction traceback (Add free gap to alignmentX)
                OptimalAlignment(Seq1, Seq2, i, j - 1, "-" + alignmentX, Seq2[j-1] + alignmentY);
            }
        }
    }
}

int main(int argc, const char *argv[]){
    ifstream fasta(argv[1]);

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
    string Seq1 = sequences[0];
    string Seq2 = sequences[1];
    // convert all characters in seq1 and seq 2 to uppercase (Case insensitive alignment program)
    transform(Seq1.begin(), Seq1.end(), Seq1.begin(), ::toupper);
    transform(Seq2.begin(), Seq2.end(), Seq2.begin(), ::toupper);

    int m = Seq1.size();
    int n = Seq2.size();

    Sc = vector(m + 1, vector<int>(n + 1));
    cout << "Alignment score:" << AlignmentScore(Seq1, Seq2, m, n) << endl;
    OptimalAlignment(Seq1, Seq2, m, n, "", "");
    return 0;
}