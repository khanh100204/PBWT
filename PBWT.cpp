#include<bits/stdc++.h>
using namespace std;    

/*      x0: 0 1 0 0 1 0
        x1: 0 1 1 0 1 1
        x2: 1 1 0 0 1 0
        x3: 0 0 1 1 1 0
        x4: 1 1 0 1 0 0
        x5: 0 0 0 1 1 0
*/

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
    vector<int> a, b;

    for (int i = 0; i < m; i++) {
        if (dk[i] > k - L) {
            if (!a.empty() && !b.empty()) {
                for (int ai : a) {
                    for (int bi : b) {
                        cout << ai << "  " << bi << "  " << k - L << "  " << k - 1<< "\n";
                    }
                }
            }
            a.clear();
            b.clear();
       
        }

        int seq = ak[i];
        if (x[seq][k] == 0)
            a.push_back(seq);
        else
            b.push_back(seq);
    }

    if (!a.empty() && !b.empty()) {
        for (int ai : a) {
            for (int bi : b) {
                cout << ai << "  " << bi << "  " << k - L << "  " << k- 1 << "\n";
            }
        }
    }
}


// Algo 4
void reportSetMaximalMatches(vector<vector<int>> X, vector<vector<int>> Ak, vector<vector<int>> Dk ) {
    int m = X.size();
    int n = X[0].size();

    for (int k = 0; k < n; k++) {
        vector<int> ak = Ak[k];
        vector<int> dk = Dk[k];

        vector<int> d(m + 2);
        d[0] = k + 1;
        for (int i = 0; i < m; i++) d[i + 1] = dk[i];
        d[m + 1] = k + 1;

        for (int i = 1; i <= m; i++) {
            int idx = ak[i - 1];
            int m = i - 1;          
            int n = i + 1;          
            bool skip = false;

            if (d[i] <= d[i + 1]) {
                while (m >= 1 && d[m + 1] >= d[i]) {    
                    int left_idx = ak[m - 1];
                    if (k < n - 1 && X[left_idx][k] == X[idx][k]) {
                        skip = true;             
                        break;
                    }
                    m--;
                }
            }
            if (skip) continue;

            if (d[i] >= d[i + 1]) {
                while (n <= m && d[n] >= d[i + 1]) {    
                    int right_idx = ak[n - 1];
                    if (k < n - 1 && X[right_idx][k] == X[idx][k]) {
                        skip = true;                    
                        break;
                    }
                    n++;
                }
            }
            if (skip) continue;

            for (int j = m + 1; j < i; j++) {
                if (k - d[i] + 1 >= L)
                    cout << idx << "  " << ak[j - 1] << "  " << d[i]<< "  " << k << "\n";
            }

            for (int j = i + 1; j < n; j++) {
                if (k - d[i] + 1 >= L)
                    cout << idx << "  " << ak[j - 1] << "  " << d[i]<< "  " << k << "\n";
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

    freopen("bin_1k_1k.txt", "r", stdin);
    freopen("bin_1k_1k.out", "w", stdout);

    vector<string> rows;
    string line;
    while (cin >> line) {
        rows.push_back(line);
    }

    
    int N = rows.size();          
    int M = rows[0].size();
    cout << "N: " << N << ", M: " << M << "\n";

    vector<vector<int>> X(M, vector<int>(N, 0));
    for (int i = 0; i < N; i++) {         
        for (int j = 0; j < M; j++) {     
            X[j][i] = rows[i][j] - '0';  
    }
}

    vector<vector<int>> Ak;
    vector<vector<int>> Dk;
    buildPBWT(X, Ak, Dk);

    // for (int i = 0; i < Ak.size(); i++) {
    //     cout << "a" << i << ": ";
    //     for (int a : Ak[i])
    //         cout << a << " ";
    //     cout << "\n";
    // }

    // cout<< "\n";
    // for (int i = 0; i < Dk.size(); i++) {
    //     cout << "d" << i << ": ";
    //     for (int d : Dk[i])
    //         cout << d << " ";
    //     cout << "\n";
    // }
    
    L = 20; // minimum match length
    for (int k = 0; k <= X[0].size(); k++) {
        reportLongMatches(X, Ak[k], Dk[k], k, L);
    }

    // reportSetMaximalMatches(X, Ak, Dk);

    return 0;
}

