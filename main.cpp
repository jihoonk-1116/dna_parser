#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cctype>
#include <algorithm>
using namespace std;
string readFastaFile(string path);
void readCsvFile(string path);
string transcribe(string dnaString);
vector<string> translate(string mrnaString);
vector<vector<string> > mycodons;

string readFastaFile(string path)
{
    string dnaString;
    fstream newfile;
    newfile.open(path, ios::in);
    if (newfile.is_open())
    {
        string line;
        getline(newfile, line);
        while (getline(newfile, line))
        {
            dnaString += line;
        }
        newfile.close();
    }
    return dnaString;
}

void readCsvFile(string path)
{
    ifstream inFile;
    string codon, Abv, AAcode, name;

    inFile.open(path);
    getline(inFile, codon, ',');
    getline(inFile, Abv, ',');
    getline(inFile, AAcode, ',');
    inFile >> name;
    while (!inFile.eof())
    {
        getline(inFile, codon, ',');
        getline(inFile, Abv, ',');
        getline(inFile, AAcode, ',');
        inFile >> name;
        vector<string> temp;
        codon.erase(std::remove_if(codon.begin(), codon.end(), ::isspace), codon.end());
        temp.push_back(codon);
        AAcode.erase(std::remove_if(AAcode.begin(), AAcode.end(), ::isspace), AAcode.end());
        temp.push_back(AAcode);
        mycodons.push_back(temp);
    }
}

string transcribe(string dna_string)
{
    string mrna_string = "";
    for (int i = 0; i < dna_string.length(); i++)
    {
        //A->U , T->A , G->C, C->G
        //ATT, ACT, ATC
        if (dna_string[i] == 'A')
            mrna_string += 'U';
        else if (dna_string[i] == 'T')
            mrna_string += 'A';
        else if (dna_string[i] == 'G')
            mrna_string += 'C';
        else if (dna_string[i] == 'C')
            mrna_string += 'G';
        else
            continue;
    }
    return mrna_string;
}

vector<string> translate(string mrna_string)
{
    //look for mRNA sequence that starts with AUG
    //UAG, UGA, UAA -> indicates the end point of a sequence
    //all sequences is ended with 'Stop'
    //AUG UGG GGA GAC ACG CGG GGC GCA CACACUCUCUCACAUAUUUUUUUAUAA
    //AUG (GAU/GAC) (GGU/GGC/GGA/GGG) (ACU/ACC/ACA/ACG) CAU/CAC UUA/UUG/CUU/CUC/CUA/CUG (AUU/AUC/AUA)
    //(UUA/UUG/CUU/CUC/CUA/CUG) AAA/AAG
    //the first : AUG GAC GGC ACU CAU UUA AUU UUA AAA UAA
    //AUG GAC GGC ACU CAU UUA AUU UUA AAA UAA
    //81
    vector<string> protein_list; //empty vector
    //UCG AAA AGU AAG ACU GAC GUU GCC CGU UAU ACA GAG ACA CAC CUA AUU UUU UUC UCA CAG
    //ACU AUC GUC GAA GAC UUG ACC A AU GGA CGG CAC UCA UUU AAU UUU AAA AUA ACU GAA UCC AGU GAU UUA UGA
    int index = 0;
    while (index <= mrna_string.length())
    {
        string result = "";
        string t = "";
        t += mrna_string[index];
        t += mrna_string[index + 1];
        t += mrna_string[index + 2];
        index = index + 1;
        //cout << t << endl;
        if (t == "AUG")
        {
            result += t;
            while (1)
            {
                t = "";
                t += mrna_string[index];
                t += mrna_string[index + 1];
                t += mrna_string[index + 2];
                //cout << t << endl;
                index = index + 3;
                //cout << index << endl;
                if (t == "UAG" || t == "UGA" || t == "UAA")
                {
                    result += t;
                    // cout << endl;
                    // cout << result << endl;
                    // cout << endl;
                    protein_list.push_back(result);
                    break;
                }
                result += t;
            }
        }
    }
    return protein_list;
}

// void print_protein_list(vector<string> list)
// {
//     for (string line : list)
//     {
//         cout << line << endl;
//     }
//     cout << list.size() << "proteins listed" << endl;
// }

int main()
{
    string dnastring = readFastaFile("ecoli.fa");
    readCsvFile("codon_table.csv");
    string mrnastring = transcribe(dnastring);
    //cout << mrnastring << endl;
    vector<string> protein_list = translate(mrnastring);
    cout << protein_list.size();
    // print_protein_list(protein_list);

    return 0;
}