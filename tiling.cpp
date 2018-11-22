
#include "tiling.h"
#include "vertex.h"

using namespace std;

class Graph
{
public:
	Vertex* s;
	Vertex* t;
	vector<vector<Vertex*>> vertexList;
	unordered_set<Vertex*> tileList;
	unordered_set<Vertex*> redSet;
	unordered_set<Vertex*> blueSet;

	Graph()
	{
		s = new Vertex();
		t = new Vertex();
	}

	void addEdge(Vertex* a, Vertex* b)
	{
		a->neighs.insert(b);
		a->weights[b] = 1;
	}
};

// Finds a (shortest according to edge length) augmenting path
// from s to t in a graph with vertex set V.
// Returns whether there is an augmenting path.
bool augmenting_path(Vertex* s, Vertex* t, 
        unordered_set<Vertex*> V, vector<Vertex*> &P)
{
        // Check that s and t aren't nullptr
        if (s == nullptr || t == nullptr)
	{
		cerr << "augmenting_path() was passed nullptr s or t." << endl;
		abort();
	}

        // Check that s and t are in the graph
        if (V.find(s) == V.end() || V.find(t) == V.end())
	{
		cerr << "augmenting_path() was passed s or t not in V." << endl;
		abort();
	}

	// Check that every vertex has valid neighs/weights.
	for (Vertex* v : V)
		for (Vertex* vn : v->neighs)
			if (v->weights.find(vn) == v->weights.end())
			{
				cerr << "augmenting_path() was passed invalid vertex." << endl;
				abort();
			}

    // Since augmenting paths should have the fewest edges,
	// not the minimum weight, run BFS.
	queue<Vertex*> Q;
	Q.push(s);

	unordered_set<Vertex*> R;
	R.clear(); 
	R.insert(s);

	unordered_map<Vertex*, Vertex*> prev;

	while (!Q.empty())
	{
		Vertex* cur = Q.front();
		Q.pop();

		for (Vertex* nei : cur->neighs)
		{
			// Must have positive edge weight
			if (cur->weights[nei] == 0)
				continue;

			if (R.find(nei) == R.end())
			{
				Q.push(nei);
				R.insert(nei);
				prev[nei] = cur; 
			}
		}
	}      

        // If BFS never reached t
        if (R.find(t) == R.end())
                return false;

        // Reconstruct shortest path backwards
        P.clear();
        P.push_back(t);
        while (P[P.size()-1] != s)
                P.push_back(prev[P[P.size()-1]]);

        // Reverse shortest path
        for (int i = 0; i < P.size()/2; ++i)
		swap(P[i], P[P.size()-1-i]);

        return true;
}

// Returns the maximum flow from s to t in a weighted graph with vertex set V.
// Assumes all edge weights are non-negative.
int max_flow(Vertex* s, Vertex* t, unordered_set<Vertex*> V)
{
	// If s or t is invalid.
        if (s == nullptr || t == nullptr)
	{
		cerr << "max_flow() was passed nullptr s or t." << endl;
		abort(); 
	}

	// If s or t is not in the vertex set.
        if (V.find(s) == V.end() || V.find(t) == V.end())
	{
		cerr << "max_flow() was passed s or t not in V." << endl;
		abort(); 
	}

	// Check that every vertex has valid neighs/weights.
	for (Vertex* v : V)
		for (Vertex* vn : v->neighs)
			if (v->weights.find(vn) == v->weights.end())
			{
				cerr << "max_flow() was passed invalid vertex." << endl;
				abort();
			}

        // Create a deep copy of V to use as the residual graph
        unordered_set<Vertex*> resV;
        unordered_map<Vertex*, Vertex*> C; // Maps vertices in V to copies in resV
        for (Vertex* vp : V)
        {
                Vertex* rp = new Vertex;
                resV.insert(rp);
                C[vp] = rp;
        }
        for (Vertex* vp : V)
                for (Vertex* np : vp->neighs)
                {
                        C[vp]->neighs.insert(C[np]);
                        C[vp]->weights[C[np]] = vp->weights[np];
                }
	// Add any missing necessary "back" edges. 
        for (Vertex* vp : V)
                for (Vertex* np : vp->neighs)
		{
			if (C[np]->neighs.find(C[vp]) == C[np]->neighs.end())
			{
				C[np]->neighs.insert(C[vp]);
				C[np]->weights[C[vp]] = 0;
			}
		}

        // Run Edmonds-Karp
        while (true)
        {
                // Find an augmenting path
                vector<Vertex*> P;
                if (!augmenting_path(C[s], C[t], resV, P))
                        break;  
                // Update residual graph
                for (int i = 0; i < P.size()-1; ++i)
                {
                        --((*(resV.find(P[i])))->weights[P[i+1]]);
                        ++((*(resV.find(P[i+1])))->weights[P[i]]);
                }
        }

        // Compute actual flow amount
        int flow = 0;
        for (Vertex* snp : C[s]->neighs)
                flow += 1 - C[s]->weights[snp];

        // Delete residual graph
        for (Vertex* vp : resV)
                delete vp;

        return flow;
}

void color(vector<vector<char>> tiles, Graph &G, unordered_map<Vertex*, string> &tileColors, int i, int j, string paintColor)
{
	Vertex* v = G.vertexList[i][j];
	Vertex* up = G.vertexList[i - 1][j];
	Vertex* down = G.vertexList[i + 1][j];
	Vertex* left = G.vertexList[i][j - 1];
	Vertex* right = G.vertexList[i][j + 1];

	// check left
	if (j != 0)
	{
		if (tiles[i][j - 1] == ' ')
		{
			tileColors[left] = paintColor;
			G.blueSet.insert(left);
			G.addEdge(v, left);
			G.tileList.insert(left);
		}
	}

	// check right
	if (j != tiles[i].size() - 1)
	{
		if (tiles[i][j + 1] == ' ')
		{
			tileColors[right] = paintColor;
			G.blueSet.insert(right);
			G.addEdge(v, right);
			G.tileList.insert(right);
		}
	}

	// check up
	if (i != 0)
	{
		if (tiles[i - 1][j] == ' ')
		{
			tileColors[up] = paintColor;
			G.blueSet.insert(up);
			G.addEdge(v, up);
			G.tileList.insert(up);
		}
	}

	// check down
	if (i != 0)
	{
		if (tiles[i + 1][j] == ' ')
		{
			tileColors[down] = paintColor;
			G.blueSet.insert(down);
			G.addEdge(v, down);
			G.tileList.insert(down);
		}
	}
}

void checkerTiles(vector<vector<char>> tiles, Graph &G, unordered_map<Vertex*, string> &tileColors)
{
	for (int i = 0; i < tiles.size(); i++)
	{
		for (int j = 0; j < tiles[i].size(); j++)
		{
		//	Vertex* v = G.vertexList[i][j];

			if (tiles[i][j] == ' ')
			{
				if (tileColors.find(G.vertexList[i][j]) == tileColors.end())
				{
					tileColors[G.vertexList[i][j]] = "red";
					G.redSet.insert(G.vertexList[i][j]);
					G.tileList.insert(G.vertexList[i][j]);
					color(tiles, G, tileColors, i, j, "blue");
				}
			}
			//delete v;
		}
	}
}

void buildMazeVector(string maze, int& numLines, Graph& G, vector<vector<char>> &arr)
{
	vector<char> row;
	vector<Vertex*> r;
	for (int i = 0; i < maze.length(); i++)
	{
		if (maze[i] != '\n')
		{
			row.push_back(maze[i]);
			r.push_back(new Vertex());
		}
		else
		{
			arr.push_back(row);
			G.vertexList.push_back(r);
			numLines++;
			row.clear();
			r.clear();
		}
	}
}

void addStart(Graph &G)
{
	for (auto it = G.redSet.begin(); it != G.redSet.end(); ++it)
		G.addEdge(G.s, *it);

	G.tileList.insert(G.s);
}

void addEnd(Graph &G)
{
	for (auto it = G.blueSet.begin(); it != G.blueSet.end(); ++it)
		G.addEdge(*it, G.t);

	G.tileList.insert(G.t);
}

bool has_tiling(string floor)
{
	Graph G;
	vector<vector<char>> tiles;
	int numLines = 0;
	unordered_map<Vertex*, string> tileColors;
	
	// Step 1: read in the string -> vector<vector<char>>
	buildMazeVector(floor, numLines, G, tiles);

	// Step 2: assign colors to each vertex (checker tiles)
	// Step 3: create the graph
	checkerTiles(tiles, G, tileColors);

	// Step 4: add start and end vertices
	addStart(G);
	addEnd(G);

	// Step 5: call max_flow
	int maxFlow = max_flow(G.s, G.t, G.tileList);

	// Step 6: return t/f
	if (maxFlow == G.redSet.size() && maxFlow == G.blueSet.size())
		return true;
	else
		return false;
}




