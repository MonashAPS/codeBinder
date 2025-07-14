
//
// Created by alexsc on 9/07/2025.
// My shitty little Aho Corasick template I copied from cp alg

#include<bits/stdc++.h>

using namespace std;


// omogen heap stuff
#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

const int K = 26; // alphabet size

/*
 * Tries are stored as lists of nodes, and edges are stored by an array of size alphabet, giving the indexes of children
 * Go and get_link functions are calculated using memoisation, calculating when necessary.
 *
 */

struct Vertex {
    int next[K]; // indexes of children in trie list
    int output = -1; //whether root -> node path matches a string
    int p = -1; //index of parent
    char pch; //char on edge from parent to this node
    int link = -1; //index of suffix link
    int go[K]; //index of transition to take. for a given i, this is next[K] if that isn't -1, else it is the suffix link
    int term = -1;
    int depth = 0;

    Vertex(int p = -1, char ch = '$', int depth = 0): p(p), pch(ch), depth(depth) {
        fill(all(next), -1);
        fill(all(go), -1);
    }
};

vector<Vertex> t(1); //trie

// add a string to the trie

void add_string(string const&s, int ind) {
    int v = 0; // index of current node in trie
    int depth = 1;
    for (char ch: s) {
        int c = ch - 'a';// get relative position of character
        if (t[v].next[c] == -1) {
            // create child node
            t[v].next[c] = t.size();
            t.emplace_back(v,ch,depth);
        }
        // move to child node
        v = t[v].next[c];

        depth++;
    }
    t[v].output = ind;
}

// at node index v, transition with character v using suffix links
int go(int v, char ch);

// get suffix link for node v;
void bfs() {
    queue<int> nodeQ;
    nodeQ.push(0);
    while (!nodeQ.empty()) {
        const int v = nodeQ.front();
        Vertex &cur = t[v];
        if (v == 0 || cur.p == 0) {
            cur.link = 0;
        }
        else {
            cur.link = go(t[cur.p].link, cur.pch);
        }
        if (t[cur.link].output != -1) {
            cur.term = cur.link;
        }
        else {
            cur.term = t[cur.link].term;
        }
        for (int child: cur.next) {
            if (child != -1) {
                nodeQ.push(child);
            }
        }

        nodeQ.pop();
    }
}


// int get_link(int v) {
//     if (t[v].link == -1) {
//         // if link doesn't exist create it
//         if (v == 0 || t[v].p == 0) {
//             // suffix link of root or direct child of root is suffix
//             t[v].link = 0;
//         }
//         else {
//             // suffix link of node is suffix link of parent, then traversing the edge from parent to this node
//             t[v].link = go(get_link(t[v].p), t[v].pch);
//         }
//     }
//     return t[v].link;
// }

int go(int v, char ch) {
    int c = ch - 'a';
    if (t[v].go[c] == -1) {
        if (t[v].next[c] != -1) {
            // take child if it exists
            t[v].go[c] = t[v].next[c];
        }
        else{
            // otherwise, traverse to suffix link and try traversing the character there
            t[v].go[c] = v == 0 ? 0: go(t[v].link, ch);
        }
    }
    return t[v].go[c];
}

int main()
{
    string txt = "abcbcbaccabcab";

    vector<string> pats  = {"a", "bc", "bca", "ac", "cc"};
    for (int i = 0; i < pats.size(); ++i) {
        add_string(pats[i],i);
    }

    bfs();

    int v = 0;
    vector<pair<int,int>> ans ;
    for (int i = 0; i < txt.size(); ++i){
        char c = txt[i];
        v = go(v,c);
        if (t[v].output != -1) {
            ans.push_back({i-t[v].depth + 1, t[v].output});
        }
        int nd = t[v].term;
        while (nd != -1) {
            ans.push_back({i-t[nd].depth +1, t[nd].output});
            nd = t[nd].term;
        }
    }

    for (auto [ind, match] : ans) {
        cout << "match for pattern: " << match << " in index : " << ind << endl;
    }
    cout << "done" << endl;
}

