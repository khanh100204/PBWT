#include<bits/stdc++.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;    

/*      x0: 0 1 0 0 1 0
        x1: 0 1 1 0 1 1
        x2: 1 1 0 0 1 0
        x3: 0 0 1 1 1 0
        x4: 1 1 0 1 0 0
        x5: 0 0 0 1 1 0
*/


// read VCF (AI generated)
vector<vector<int>> readVCF(const string& filename) {

    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Cannot open VCF file: " << filename << endl;
        exit(1);
    }

    vector<vector<int>> X;

    string line;
    int numSamples = 0;

    while (getline(file, line)) {

        if (line.rfind("##", 0) == 0)
            continue;

        if (line.rfind("#CHROM", 0) == 0) {

            stringstream ss(line);
            string field;
            vector<string> fields;

            while (ss >> field)
                fields.push_back(field);

            numSamples = fields.size() - 9;

            X.resize(numSamples * 2);

            continue;
        }

        if (line.empty())
            continue;

        stringstream ss(line);
        string field;
        for (int i = 0; i < 9; i++)
            ss >> field;

        for (int sample = 0; sample < numSamples; sample++) {

            string genotype;
            ss >> genotype;

            // In a normal VCF this could be:
            // 0|1
            // or 0|1:35:99:...
            // We only care about GT.
            size_t colon = genotype.find(':');
            if (colon != string::npos)
                genotype = genotype.substr(0, colon);

            size_t bar = genotype.find('|');

            if (bar == string::npos) {
                cerr << "Expected phased genotype, got: "
                     << genotype << endl;
                exit(1);
            }

            int allele1 = genotype[0] - '0';
            int allele2 = genotype[bar + 1] - '0';

            X[2 * sample].push_back(allele1);
            X[2 * sample + 1].push_back(allele2);
        }
    }

    return X;
}


int L; // minimum match length

// Algo 1 and 2
void buildPBWT(const vector<vector<int>>& x, vector<vector<int>>& Ak, vector<vector<int>>& Dk) {
    int m = x.size();
    int n = x[0].size();

    vector<int> yk(m);
    vector<int> dk(m, 0);
    for (int i = 0; i < m; i++) {
        yk[i] = i;
    }
    Ak.push_back(yk);
    Dk.push_back(dk);
    for (int k = 0 ; k < n ; k++) {
        int p = k + 1, q = k + 1;
        vector<int> a, b, d, e;

        
/*      x0: 0 1 0 0 1 0
        x1: 0 1 1 0 1 1
        x2: 1 1 0 0 1 0
        x3: 0 0 1 1 1 0
        x4: 1 1 0 1 0 0
        x5: 0 0 0 1 1 0
*/

        for (int i = 0 ; i < m ; i++) {
            int real_index = yk[i];
            if (dk[i] > p) p = dk[i];
            if (dk[i] > q) q = dk[i];

            
            if (x[real_index][k] == 0) {
                a.push_back(real_index);
                d.push_back(p);
                p = 0;  
            } 
            else {
                b.push_back(real_index);
                e.push_back(q);
                q = 0;
            }
        }
        yk.clear();
        yk.insert(yk.end(), a.begin(), a.end());
        yk.insert(yk.end(), b.begin(), b.end());
        Ak.push_back(yk);
   
        dk.clear();
        dk.insert(dk.end(), d.begin(), d.end());
        dk.insert(dk.end(), e.begin(), e.end());
        Dk.push_back(dk); 
   }

    return;
}

// Algo 3
void reportLongMatches(vector<vector<int>> x, vector<int> ak, vector<int> dk, int k, int L ) {
    int m = x.size();
    vector<pair<int,int>> a, b;

    auto reportBlock = [&]() {
        if (a.empty() || b.empty()) return;
        for (auto [ai, posA] : a) {
            for (auto [bi, posB] : b) {
                int left = min(posA, posB);
                int right = max(posA, posB);
                int start = 0;

                for (int j = left + 1; j <= right; j++) {
                    start = max(start, dk[j]);
                }
                // report
                cout << min(ai,bi) << "  " << max(ai,bi) << "  " << start << "  " << k - 1<< "\n";

            }
        }
    };


    for (int i = 0; i < m; i++) {
        if (dk[i] > k - L) {
            reportBlock();
            a.clear();
            b.clear();
       
        }

        int seq = ak[i];
        if (x[seq][k] == 0)
            a.push_back({seq,i});
        else
            b.push_back({seq,i});
    }
    // last flush
    reportBlock();
}


// Algo 4
void reportSetMaximalMatches(vector<vector<int>> X, vector<vector<int>> Ak, vector<vector<int>> Dk ) {
    int M = X.size();          
    int N = X[0].size();       

    for (int k = 0; k <= N; k++) {

        const vector<int>& ak = Ak[k];

        vector<int> d = Dk[k];

        d[0] = k + 1;
        d.push_back(k + 1);

        bool forcedEnd = (k == N);

        for (int i = 0; i < M; i++) {

            int mh = ak[i];

            int left = i - 1;
            int right = i + 1;

            bool skip = false;

            if (d[i] <= d[i + 1]) {

                while (d[left + 1] <= d[i]) {

                    int neighbor = ak[left];

                    
                    if (!forcedEnd && X[neighbor][k] == X[mh][k]) {
                        skip = true;
                        break;
                    }

                    left--;
                }
            }

            if (skip) {
                continue;
            }

            if (d[i] >= d[i + 1]) {

                while (d[right] <= d[i + 1]) {

                    int neighbor = ak[right];

                    if (!forcedEnd && X[neighbor][k] == X[mh][k]) {
                        skip = true;
                        break;
                    }

                    right++;
                }
            }

            if (skip) {
                continue;
            }

            int startLeft = d[i];
            int lengthLeft = k - startLeft;

            if (lengthLeft > 0) {

                for (int j = left + 1; j < i; j++) {

                    int neighbor = ak[j];

            
                    cout << mh << "  " << neighbor << "  " << k + 1 << "  " << lengthLeft << "\n";
                }
            }

           
            int startRight = d[i + 1];
            int lengthRight = k - startRight;
 
            if (lengthRight > 0) {

                for (int j = i + 1; j < right; j++) {

                    int neighbor = ak[j];

                    cout << mh << "  "  << neighbor << "  " << k + 1 << "  " << lengthRight << "\n";
                }
            }
        }
    }
}

// vector<vector<int>> X = {
    //         {0, 1, 0, 0, 1, 1},
    //         {0, 1, 1, 0, 1, 1},
    //         {1, 1, 0, 0, 1, 0},
    //         {0, 0, 1, 1, 1, 0},
    //         {1, 1, 0, 1, 0, 0},
    //         {0, 0, 0, 1, 1, 0}  
    //     };

int main(){

    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    // freopen("bin_1k_1k.txt", "r", stdin);
    freopen("algo4.out", "w", stdout);

    vector<vector<int>> X = readVCF("Algo4/s100i50.rnd.vcf");

    // /* easy test */
    // vector<vector<int>> X = {
    //         {0, 1, 0, 0, 1, 1},
    //         {0, 1, 1, 0, 1, 1},
    //         {1, 1, 0, 0, 1, 0},
    //         {0, 0, 1, 1, 1, 0},
    //         {1, 1, 0, 1, 0, 0},
    //         {0, 0, 0, 1, 1, 0}  
    //     };


    vector<vector<int>> Ak;
    vector<vector<int>> Dk;
    buildPBWT(X, Ak, Dk);


    
    L = 20; // minimum match length
    // for (int k = 0; k <= X[0].size(); k++) {
    //     reportLongMatches(X, Ak[k], Dk[k], k, L);
    // }

    reportSetMaximalMatches(X, Ak, Dk);

    return 0;
}

