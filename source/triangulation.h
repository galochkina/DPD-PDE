#ifndef TRIANGULATION_H
#define	TRIANGULATION_H

typedef Point VecF2;

//*******************************************************************************************************************************************************
class Segment
{
public:
	VecF2 v0, v1;		// two endpoints
	int is_inserted;	// inserted in trapezoidation yet?
	int root0, root1;	// root nodes in Q
	int next;			// next logical segment
	int prev;			// previous segment

private:

};

//*******************************************************************************************************************************************************
class Trap //trapezoid
{
public:
	int lseg, rseg;		// two adjoining segments
	Point hi, lo;		// max/min y-values
	int u0, u1;
	int d0, d1;
	int sink;			// pointer to corresponding in Q
	int usave, uside;	// I forgot what this means
	int state;

private:

};

//*******************************************************************************************************************************************************
class Node //node in the query structure
{
public:
	int nodetype;		// Y-node or S-node
	int segnum;
	Point yval;
	int trnum;
	int parent;			// doubly linked DAG
	int left, right;	// children

private:

};

//*******************************************************************************************************************************************************
class Monchain
{
public:
	int vnum;
	int next;			// circularly linked list
	int prev;			// describing the monotone
	int marked;			// polygon

private:

};

//*******************************************************************************************************************************************************
class Vertexchain
{
public:
	Point pt;
	int vnext[4];		// next vertices for the 4 chains
	int vpos[4];		// pos of v in the 4 chains
	int nextfree;

private:

};

//*******************************************************************************************************************************************************
class Triang
{
public:
	vector<Node> qs;
	vector<Trap> tr;
	vector<Segment> seg;

	Triang(){}

	void execute(int ncontours)
	{
		int npoints;
		int ccount = 0;
		while(ccount < ncontours)
		{
			npoints = cntr[ccount];
		}
	}
	

private:


};

#endif