#include<bits/stdc++.h>
using namespace std;    

/*      x0: 0 1 0 0 1 0
        x1: 0 1 1 0 1 1
        x2: 1 1 0 0 1 0
        x3: 0 0 1 1 1 0
        x4: 1 1 0 1 0 0
        x5: 0 0 0 1 1 0
*/

struct Node {
    int ID;                 // ak[i]
    int d;                  // dk[i]

    Node* above;
    Node* below;

    Node* u;                // follow 0 
    Node* v;                // follow 1 

    bool sentinel = false;

    Node* w(int bit) {// w[0] = u, w[1] = v
        return (bit == 0 ? u : v);
    }

    void setW(int bit, Node* p) {
        if (bit == 0) u = p;
        else v = p;
    }
};

class DPBWT {
public:

    vector<vector<int>> X;

    int N;                      // number of sites

    // one linked list for each a0 ... aN
    vector<Node*> head;
    vector<Node*> tail;

    // byID[id][k] = node representing xid in column k
    vector<vector<Node*>> byID;
    DPBWT(const vector<vector<int>>& X0, const vector<vector<int>>& Ak, const vector<vector<int>>& Dk) {    
        X = X0;
        N = X[0].size();

        int M = X.size();

        head.resize(N + 1);
        tail.resize(N + 1);

        byID.assign(M, vector<Node*>(N + 1));

        for (int k = 0; k <= N; k++) {

            head[k] = new Node();
            tail[k] = new Node();

            head[k]->sentinel = true;
            tail[k]->sentinel = true;

            head[k]->below = tail[k];
            tail[k]->above = head[k];

            Node* last = head[k];

            for (int i = 0; i < M; i++) {

                Node* n = new Node();

                n->ID = Ak[k][i];
                n->d  = Dk[k][i];

                n->above = last;
                n->below = tail[k];

                last->below = n;
                tail[k]->above = n;

                last = n;

                byID[n->ID][k] = n;
            }
        }

        buildUV();
    };

    void buildUV() {
        for (int k = 0; k < N; k++) {

            //Find boundary
            Node* firstOne = tail[k + 1];

            for (Node* p = head[k + 1]->below; p != tail[k + 1]; p = p->below) {

                if (X[p->ID][k] == 1) {
                    firstOne = p;
                    break;
                }
            }

            /*  end:            
                follow 0 -> end of zero block??
                follow 1 -> bottom??
            */

            tail[k]->u = firstOne;
            tail[k]->v = tail[k + 1];

            Node* next0 = firstOne;
            Node* next1 = tail[k + 1];

            /*
                next0 = first 0 at/below current node
                next1 = first 1 at/below current node
            */

            for (Node* p = tail[k]->above; p != head[k]; p = p->above) {

                int bit = X[p->ID][k];

                if (bit == 0)
                    next0 = byID[p->ID][k + 1];
                else
                    next1 = byID[p->ID][k + 1];

                p->u = next0;
                p->v = next1;
            }
        }
    }

    void insertBefore(Node* pos, Node* n) {
        Node* above = pos->above;

        n->above = above;
        n->below = pos;

        above->below = n;
        pos->above = n;
    }


    void insertHaplotype(const vector<int>& z) {

            int newID = X.size();

            X.push_back(z);

            byID.push_back(
                vector<Node*>(N + 1, nullptr)
            );

            vector<Node*> zn(N + 1);
            vector<Node*> t(N + 1);


            // --------------------------------
            // FORWARD SWEEP
            // --------------------------------

            zn[0] = new Node();
            zn[0]->ID = newID;

            insertBefore(t[0], zn[0]);

            byID[newID][0] = zn[0];


            for (int k = 0; k < N; k++) {

                Node* tk = t[k];

                Node* nextT =
                    tk->w(z[k]);

                Node* oppositeTarget =
                    tk->w(1 - z[k]);


                zn[k + 1] = new Node();
                zn[k + 1]->ID = newID;

                insertBefore(nextT, zn[k + 1]);

                byID[newID][k + 1] = zn[k + 1];

                t[k + 1] = nextT;
                zn[k]->setW(
                    z[k],
                    zn[k + 1]
                );
                zn[k]->setW(
                    1 - z[k],
                    oppositeTarget
                );

                Node* p = zn[k]->above;

                while (p != head[k] &&
                    X[p->ID][k] != z[k]) {

                    p->setW(
                        z[k],
                        zn[k + 1]
                    );

                    p = p->above;
                } 
        }
        // BACKWARD SWEEP: divergence
        // --------------------------------

        int upStart = N;
        int downStart = N;

        for (int k = N; k >= 0; k--) {

            upStart = min(upStart, k);
            downStart = min(downStart, k);

            Node* up =
                zn[k]->above;

            Node* down =
                zn[k]->below;


            // divergence between ABOVE and z

            if (up == head[k]) {

                zn[k]->d = k;

            } else {

                while (
                    upStart > 0 &&
                    z[upStart - 1]
                        ==
                    X[up->ID][upStart - 1]
                ) {
                    upStart--;
                }

                zn[k]->d = upStart;
            }


            // divergence between z and BELOW

            if (down != tail[k]) {

                while (
                    downStart > 0 &&
                    z[downStart - 1]
                        ==
                    X[down->ID][downStart - 1]
                ) {
                    downStart--;
                }

                down->d = downStart;
            }
        }
    }

    void longMatchQueryAlg3(const vector<int>& z, int L) {

        if (L <= 0 || L > N)
            return;


        int M = X.size();

        // SWEEP 1
        // Virtually insert z.

        vector<Node*> t(N + 1);

        t[0] = tail[0];

        for (int k = 0; k < N; k++) {

            t[k + 1] =
                t[k]->w(z[k]);
        }

        // SWEEP 2
        // virtual divergence:
        // zd[k] = divergence(z, sequence above z)
        // bd[k] = divergence(z, sequence below z)

        vector<int> zd(N + 1);
        vector<int> bd(N + 1);

        int zStart = N;
        int bStart = N;


        for (int k = N; k >= 0; k--) {

            zStart = min(zStart, k);
            bStart = min(bStart, k);

            Node* below = t[k];
            Node* above = below->above;


            if (above != head[k]) {

                while (
                    zStart > 0 &&
                    z[zStart - 1]
                        ==
                    X[above->ID][zStart - 1]
                ) {
                    zStart--;
                }

                zd[k] = zStart;

            } else {

                zd[k] = k;
            }


            if (below != tail[k]) {

                while (
                    bStart > 0 &&
                    z[bStart - 1]
                        ==
                    X[below->ID][bStart - 1]
                ) {
                    bStart--;
                }

                bd[k] = bStart;

            } else {

                bd[k] = k;
            }
        }


        // SWEEP 3
        // [f,g) = block of sequences currently
        //          matching z for >= L


        vector<int> matchStart(M, 0);

        Node* f = t[L - 1];
        Node* g = t[L - 1];


        auto reportRange =
            [&](Node* first,
                Node* last,
                int endExclusive) {

            for (
                Node* p = first;
                p != last;
                p = p->below
            ) {
                cout << "z vs x" << p->ID << " : [" << matchStart[p->ID]  << ", " << endExclusive << ")\n";
            }
        };


        for (int k = L - 1;
            k < N;
            k++) {

            int bit = z[k];


            // which existing matches END here?
            // they take the opposite bit???

            Node* fEnd =
                f->w(1 - bit);

            Node* gEnd =
                g->w(1 - bit);

            f = f->w(bit);
            g = g->w(bit);

            // Report [start, k)

            reportRange(
                fEnd,
                gEnd,
                k
            );


            int col =
                k + 1;

            int threshold =
                col - L;

            if (f == g) {

                if (
                    f->above != head[col] &&
                    zd[col] == threshold
                ) {

                    f = f->above;

                    matchStart[f->ID] =
                        threshold;
                }


                if (
                    g != tail[col] &&
                    bd[col] == threshold
                ) {

                    matchStart[g->ID] =
                        threshold;

                    g = g->below;
                }
            }

            if (f != g) {


                while (
                    f->above != head[col] &&
                    f->d <= threshold
                ) {

                    f = f->above;

                    matchStart[f->ID] =
                        threshold;
                }

                /*expand downward*/

                while (
                    g != tail[col] &&
                    g->d <= threshold
                ) {

                    matchStart[g->ID] =
                        threshold;

                    g = g->below;
                }
            }
        }


        reportRange(
            f,
            g,
            N
        );
    }

};



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

int main() {

    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    vector<vector<int>> X = {
        {0, 1, 0, 0, 1, 1},
        {0, 1, 1, 0, 1, 1},
        {1, 1, 0, 0, 1, 0},
        {0, 0, 1, 1, 1, 0},
        {1, 1, 0, 1, 0, 0},
        {0, 0, 0, 1, 1, 0}
    };

    vector<vector<int>> Ak;
    vector<vector<int>> Dk;

    buildPBWT(X, Ak, Dk);

    for (int i = 0; i < Ak.size(); i++) {
        cout << "a" << i << ": ";

        for (int a : Ak[i])
            cout << a << " ";

        cout << "\n";
    }

    cout << "\n";

    for (int i = 0; i < Dk.size(); i++) {
        cout << "d" << i << ": ";

        for (int d : Dk[i])
            cout << d << " ";

        cout << "\n";
    }


    DPBWT dpbwt(X, Ak, Dk);

    vector<int> z = {
        0, 1, 0, 1, 1, 0
    };

    int L = 2;

    cout << "\nLong matches with z:\n";

    dpbwt.longMatchQueryAlg3(z, L);


    return 0;
}