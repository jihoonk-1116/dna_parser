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
    vector<string> protein_list; //empty vector
    int index = 0;               //string index
    while (index <= mrna_string.length())
    {
        //cout << mrna_string << endl;
        int aug_pos = mrna_string.find("AUG"); //return aug index
        index = aug_pos + 3;                   //after aug
        string result = "M";
        while (1)
        {
            string t = "";
            t += mrna_string[index];
            t += mrna_string[index + 1];
            t += mrna_string[index + 2];
            index = index + 3;
            if (index > mrna_string.length())
                break;
            for (int i = 0; i < mycodons.size(); i++)
            {
                if (mycodons[i][0] == t)
                {
                    result += mycodons[i][1];
                }
            }
            if (t == "UAG" || t == "UGA" || t == "UAA")
            {
                protein_list.push_back(result);
                mrna_string = mrna_string.substr(index);
                break;
            }
        }
    }
    return protein_list;
}

int main()
{
    string dnastring = readFastaFile("ecoli.fa");
    readCsvFile("codon_table.csv");
    string mrnastring = transcribe(dnastring);
    vector<string> protein_list = translate(mrnastring);
    cout << protein_list.size();
    // print_protein_list(protein_list);

    return 0;
}