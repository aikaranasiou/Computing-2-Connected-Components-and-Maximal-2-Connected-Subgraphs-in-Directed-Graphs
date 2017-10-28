#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <string.h>
#include "rfw_timer.h"
#include <math.h>

#include "utilities.h"

//using namespace std;

#define MAXLINE       100   /* max length of input line */
char line[MAXLINE]; /* stores input line */

int n, m; /* number of nodes, arcs */
int *input_file; // stores input arcs
FILE *output;

// dynamic list
int *First;
int *Target;
int *Next;
int *Before;
int *r_First;
int *r_Target;
int *r_Next;
int *r_Before;
int *marked_S1;
int *marked_S2;
int current_pos;
int r_current_pos;
int starting_vertex;

int *label2pre;
int *pre2label;

int *parent;
int *label2pre_2;
bool *marked_Cond;
bool *in_small_comp;

int *w; // weight of node i 

int next,next_2;
int delta;
int rand_node;
int F_one, F_two;
int check; // check=2*delta+1 

int s, p_S1,p_S2,p,  number_of_nodes_of_1egeOut,idf;
int number_of_components;
int *stack, stack_top; 
bool *stack_Cond;
bool *marked_S1_Cond;
bool *marked_S2_Cond;

///// SCC


int maxRecursionDepth, max2ECSSize, min2ECSSize, sum2ECSSize, numOf2ECS, numOfRecursionCalls;
int max1EOutSize, min1EOutSize, sum1EOutSize, numOf1EOut;


int *SCCs, *SCCs_first, *SCC_ID, SCC_next_pos; // structures to store strongly connected components

bool* SCC_isInStack;

int *temp_stack, temp_stack_bottom, temp_stack_top;
bool *temp_inserted;


bool is_export_output = true;
string algorithm_code = "2E-CHILP";
float avg2eccSize = 0;
string inputFileName = "";
string preText = "";
string testDate = "";
double takenTime = 0;

void ShowOutput();
void ExportOutput();

inline void insert_arc_to_graph(int x,int y){
    if(First[x]==0){
        Target[current_pos]=y;
        First[x]=current_pos++;
    }
    else{
        Target[current_pos]=y;
        Next[current_pos]=First[x];
	Before[First[x]]=current_pos;
        First[x]=current_pos++;
    }
    if(r_First[y]==0){
        r_Target[r_current_pos]=x;
        r_First[y]=r_current_pos++;
    }
    else{
        r_Target[r_current_pos]=x;
        r_Next[r_current_pos]=r_First[y];
	r_Before[r_First[y]]=r_current_pos;
        r_First[y]=r_current_pos++;
    }
}

inline void delete_arc_from_graph(int x,int y, int pos){
    if(First[x]==pos){
        First[x]=Next[pos];
	Before[Next[pos]] = 0;
    }
    else if(Next[pos]==0) {
        Next[Before[pos]] = Next[pos];
    }
    else{
        Next[Before[pos]] = Next[pos];
	Before[Next[pos]] = Before[pos];
    }
    
    if(r_First[y]==pos){
        r_First[y]=r_Next[pos];
	r_Before[r_Next[pos]] = 0;
    }
    else if(r_Next[pos]==0) {
        r_Next[r_Before[pos]] = r_Next[pos];
    }
    else{
        r_Next[r_Before[pos]] = r_Next[pos];
	r_Before[r_Next[pos]] = r_Before[pos];
    }
}

/* Put all reachable vertices into the temporary stack*/
void fillReachableStack(int source, int* stack){
    temp_stack_top = temp_stack_bottom = 0;
    int temp_vertex;

    stack[temp_stack_top++] = source;
    temp_inserted[source] = true;

    while(temp_stack_top!= temp_stack_bottom){
        temp_vertex = stack[temp_stack_bottom++];
	
	
	int j = First[temp_vertex];
	while(j!=0){
	      int v = Target[j];
	      if(!temp_inserted[v]){
		  stack[temp_stack_top++] = v;
		  temp_inserted[v] = true;
	      }
	      j=Next[j];
	}
	
	j = r_First[temp_vertex];
	while(j!=0){
	      int v = r_Target[j];
	      if(!temp_inserted[v]){
		  stack[temp_stack_top++] = v;
		  temp_inserted[v] = true;
	      }
	      j=r_Next[j];
	}
    }
    for(int i = 0; i < temp_stack_top; i++) temp_inserted[stack[i]] = false;
}

// initializes vector only for vertices currently in the graph
void initReachableToFalse(bool *vector){
    for(int i = 0; i < temp_stack_top; i++){
        vector[temp_stack[i]] = false;
    }
}

// initializes vector only for vertices currently in the graph
void initReachableToInt(int *vector, int value){
    for(int i = 0; i < temp_stack_top; i++){
        vector[temp_stack[i]] = value;
    }
}

//////// FIND STRONG BRIDGES ////////


int *LowLink;

// computes thw strongly connected components of a graph
void StrongConnect(int k, int *N, int *currentSccID, int *SCC_next_pos, bool* isInStack){

    int w;
    LowLink[k] = label2pre[k] = ++(*N);
    stack[++stack_top] = k;
    isInStack[k] = true;

    int j = First[k];
    while(j!=0){
	  int v = Target[j];
	  
	  if(!label2pre[v]){ // tree arc
	      StrongConnect(v, N, currentSccID, SCC_next_pos, isInStack);
	      LowLink[k] = min(LowLink[k], LowLink[v]);
	      
	  }
	  else if(label2pre[v] < label2pre[k]){ // front or cross arc
	      if(isInStack[v]) LowLink[k] = min(LowLink[k], label2pre[v]);
	  }
	  j=Next[j];
    }
    if(LowLink[k] == label2pre[k]){ // k is the root of a component
        while( (stack_top > 0) && (label2pre[w = stack[stack_top]] >= label2pre[k]) ){
            SCC_ID[w] = (*currentSccID);
            isInStack[w] = false;
            SCCs[(*SCC_next_pos)++] = w;
            stack_top--;
        }
        SCCs_first[++(*currentSccID)] = (*SCC_next_pos);
    }
}

// initializes structures and maintains the SCC that contains source or the SCCs containing the neighbors of source
void connectedComponents(int *currentSccID){
    int i;

    bool *isInStack = temp_inserted;

    // find vertices that are in the same SCC with v 
    initReachableToInt(label2pre,0);
    initReachableToInt(SCC_ID,0);

    int N = 0;
    int SCC_next_pos = 1;

    // mark the beginning of vertices of the first SCC
    SCCs_first[*currentSccID] = SCC_next_pos;

    for(i = 0; i < temp_stack_top; i++){
        if( !label2pre[temp_stack[i]] ){
            // compute the SCC containing source
            StrongConnect(temp_stack[i], &N, currentSccID, &SCC_next_pos, isInStack);
        }
    }

    initReachableToFalse(temp_inserted);

/*    // print the SCCs computed
    printf("Number of strongly connected components : %d\n", (*currentSccID)-1);
    i = 1;
    while(i < (*currentSccID)){
        fprintf(stdout,"In SCC with ID %d (%d vertices):",i,SCCs_first[i+1]-SCCs_first[i]);
        for(int j = SCCs_first[i]; j < SCCs_first[i+1]; j++){
            fprintf(stdout," %d",SCCs[j]);
        }
        fprintf(stdout,"\n");
        i++;
    }
*/
}


/* depth-first search from vertex k
   stores parents and labels */
int *dfs_parent;

void DFS(int k, int *N, int *parent,bool *skip) {
  
    label2pre[k] = ++(*N);
    pre2label[*N] = k;

    int j = First[k];
    while(j!=0){
	int v = Target[j];
	if ( !label2pre[v] ){
            dfs_parent[v] = -k; // use negative for vertices that have not been processed
            DFS(v, N, parent, skip);
            parent[label2pre[v]] = label2pre[k];
        }
        else if ( dfs_parent[v] > 0 ) skip[v] = true; // (k,v) is forward or cross arc
        
	j=Next[j];
    }

    dfs_parent[k] = -dfs_parent[k]; // k is processed
}

inline void rcompress(int v, int *parent, int *semi, int *label, int c) {
    int p;
    if ((p = parent[v]) > c) {
        rcompress(p, parent, semi, label, c);
        if (semi[label[p]] < semi[label[v]]) label[v] = label[p];
        parent[v] = parent[p];
    }
}


bool reverse;

/* depth-first search of dominator tree from node k; assigns the root of each subtree */
void SubTreeDFS(int k, int root, int *d_adj, int *d_adj_first, int *SubTree, bool *marked){
    int new_root = root;
    if(marked[k]) new_root = k;

    SubTree[pre2label[k]] = new_root;
    //if(SubTree[pre2label[k]]!= 1 && reverse == false) printf("Vertex %d is dominated by a strong bridge in G\n",pre2label[k]);
    //if(SubTree[pre2label[k]]!= 1 && reverse == true) printf("Vertex %d is dominated by a strong bridge in G^R\n",pre2label[k]);

    for (int j = d_adj_first[k]; j < d_adj_first[k+1]; j++) {
        int v = d_adj[j];
        SubTreeDFS(v, new_root, d_adj, d_adj_first, SubTree, marked);
    }
}

int *dlabel2pre, *dsize;

/* depth-first search of dominator tree from node k; assigns preorder numbers and sizes */
void DomDFS(int k, int *dN, int *d_adj, int *d_adj_first,int *dpre2label)
{
    dlabel2pre[pre2label[k]]  = ++(*dN);
    dpre2label[*dN] = pre2label[k];

    for (int j=d_adj_first[k]; j<d_adj_first[k+1]; j++) {
        int v = d_adj[j];
        DomDFS(v, dN, d_adj, d_adj_first, dpre2label);
        dsize[pre2label[k]] += dsize[pre2label[v]];
    }
}

/* constructs the dominator tree as a graph and perfors a DFS*/
inline void ConstructDomTree(int verticesNum, int *Dout, int *DfirstOut, int *dom,int *dpre2label) {
    int bsize = verticesNum + 2;
    int *DnextOut  = new int [bsize];

    for (int i = 0; i <= verticesNum; i++) {
        DfirstOut[i] = 0;
        dpre2label[i] = 0;
    }

    DfirstOut[verticesNum+1] = DnextOut[verticesNum+1] = 0;
    initReachableToInt(dsize, 1);
    initReachableToInt(dlabel2pre, 0);

    int v;
    for (v = 1; v <= verticesNum; v++) {
        if (!dom[v] || v == dom[v]) continue; // skip unreachable vertices
        DfirstOut[dom[v]+1]++;
    }

    for (int i = 1; i <= verticesNum+1; i++) {
        DfirstOut[i] += DfirstOut[i-1];
        DnextOut[i] = DfirstOut[i];
    }

    for (v = 1; v <= verticesNum; v++) {
        if (!dom[v] || v == dom[v]) continue; // skip unreachable vertices
        Dout[DnextOut[dom[v]]++] = v;
    }

    delete [] DnextOut;
    //printGraph(Dout, DfirstOut, n);
}


int max_bridge_size;
int max_bridge;

// finds all the strong bridges of a given graph (dominatos required)
inline int findBridges(int verticesNum, int* dom, int *SubTree, bool *skip) {
    int k;
    int skipped = 0, num_of_bridges = 0;
    bool foundBridge;

    int *dpre2label = new int[verticesNum+1];
    int *Dout      = new int [verticesNum+2];
    int *DfirstOut = new int [verticesNum+2];
    bool *marked = new bool [verticesNum+1];

    // construct dominator tree
    ConstructDomTree(verticesNum, Dout, DfirstOut, dom, dpre2label);

    int dN = 0;
    DomDFS(1, &dN, Dout, DfirstOut, dpre2label);

    for(int i = 0; i <=verticesNum; i++) marked[i] = false;
    
    max_bridge_size =0;
    max_bridge =0;

    for(int i = verticesNum; i>1; i--) // process in reverse preorder
    {
        k = pre2label[i]; // process ith vertex in preorder
        if (skip[k] ) {
            //printf("vertex %d skipped\n",k);
            skipped++;
            continue;
        }
        
        //printf("process node %d (pre = %d), idom %d, parent %d\n",k,i,pre2label[dom[i]],dfs_parent[k]);

        foundBridge = true;


	int j = r_First[k];
	while(j!=0){
	    int v = r_Target[j];
	    //printf("\t arc (%d,%d)\n",u,k);
	    if ( v == dfs_parent[k] ){
	      j=r_Next[j];
	      continue;
	    }

	    // check if k does not dominates u
	    if ( ( dlabel2pre[v] < dlabel2pre[k] ) || ( dlabel2pre[v] >= dlabel2pre[k] + dsize[k] ) ) {
		foundBridge = false;
		break;
	    }
	    j=r_Next[j];
	}

        //printf("done node %d, bridge %d\n",k,foundBridge);
        if (foundBridge) {
	    //if(reverse)	printf("edge (%d,%d) is a bridge in G\n", k,dfs_parent[k]);
	    //else 	printf("edge (%d,%d) is a bridge in G^R\n", dfs_parent[k],k);
	    if(dsize[k]>max_bridge_size && (verticesNum-dsize[k]) > max_bridge_size){
		  if(dsize[k] < (verticesNum-dsize[k]))	max_bridge_size = dsize[k];
		  else 					max_bridge_size = verticesNum-dsize[k];
		  max_bridge = k;
	    }
            num_of_bridges++;
            //marked[label2pre[k]] = true;
        }
    }
    if(max_bridge_size > 0) marked[label2pre[max_bridge]] = true;

    SubTreeDFS(1, 1, Dout, DfirstOut, SubTree, marked);

    delete [] dpre2label;
    delete [] Dout;
    delete [] DfirstOut;
    delete [] marked;

    return num_of_bridges;
}

/*computed dominators of a flow graph rooted at r1*/
int slt(int r1, int num_of_vertices, int *SubTree) {

    int bsize = num_of_vertices + 1;
    int *parent = new int [bsize];
    int *semi = new int [bsize];
    int *label= new int [bsize];
    int *ubucket = new int [bsize];
    int *dom = new int [bsize];
    int i;

    dfs_parent = marked_S1;
    dlabel2pre = marked_S2;
    dsize = label2pre_2;
    
    for (i = num_of_vertices; i >= 0; i--) {
        label[i] = semi[i] = i;
        ubucket[i] = 0;
        dom[i] = 0;
    }

    // initialize matrices only for the reachable vertices
    //initReachableToInt(dfs_parent, 0);
    initReachableToInt(label2pre, 0);
    
    
    bool *skip = marked_Cond;
    initReachableToFalse(skip);

    int N = 0;

    parent[1] = 1;
    DFS(r1, &N, parent, skip);

    // process the vertices in reverse preorder
    for (i = N; i > 1; i--) {
        /*---------------------
         | process i-th bucket
         *--------------------*/
        for (int v = ubucket[i]; v; v = ubucket[v]) {
            rcompress(v, parent, semi, label, i);
            int u = label[v]; // u is a vertex with min sdom in (sdom[v],v]
            dom[v] = (semi[u] < semi[v]) ? u : i;
        }
        //no need to empty the bucket

        /*---------------------------------------------
         | check incoming arcs, update semi-dominators
         *--------------------------------------------*/
        int k = pre2label[i];
		
	int j = r_First[k];
	while(j!=0){
	    int v = label2pre[r_Target[j]];
            if (v) {
                int u;
                if (v <= i) {
                    u = v;
                }//v is an ancestor of i
                else {
                    rcompress(v, parent, semi, label, i);
                    u = label[v];
                }

                if (semi[u] < semi[i]) {
                    semi[i] = semi[u];
                }
            }
	    j=r_Next[j];
	}
	
        /*---------------------------
         | process candidate semidom
         *--------------------------*/
        int s = semi[i];
        if (s != parent[i]) { //if semidominator n not parent: add i to s's bucket
            ubucket[i] = ubucket[s];
            ubucket[s] = i;
        } else {
            dom[i] = s; //semidominator is parent: s is a candidate dominator
        }
    }
    
    /*------------------
     | process bucket 1
     *-----------------*/
    for (int v = ubucket[1]; v; v = ubucket[v]) dom[v] = 1;
    
    /*---------------
     | recover idoms
     *--------------*/

    dom[1] = 1;
    for ( i = 2; i <= N; i++) {
        if (dom[i] != semi[i]) dom[i] = dom[dom[i]]; //make relative absolute // stores preorder
    }
    
    // find bridges ...
    int num_of_bridges = findBridges(N, dom, SubTree, skip);
    
    delete [] parent;
    delete [] semi;
    delete [] label;
    delete [] ubucket;
    delete [] dom;

    return num_of_bridges;
}

int *SubTreeForward,*SubTreeReverse;

int deleteBridges(int v, int num_of_vertices){
    
    // find vertices that are in the same SCC with v 
    fillReachableStack(v, temp_stack);

    SubTreeForward = w;
    SubTreeReverse = LowLink;
    
    initReachableToInt(SubTreeForward,0);
    initReachableToInt(SubTreeReverse,0);

    reverse = false;
    int BridgesNum = slt(v, num_of_vertices, SubTreeForward);

    int best_edge_disconnects_F = max_bridge_size;
    
    // swap adjacency lists
    int *temp_adj = First;
    First = r_First;
    r_First = temp_adj;
    
    temp_adj = Target;
    Target = r_Target;
    r_Target = temp_adj;
    
    temp_adj = Next;
    Next = r_Next;
    r_Next = temp_adj;
    
    temp_adj = Before;
    Before = r_Before;
    r_Before = temp_adj;
  
    reverse = true;
    BridgesNum += slt(v, num_of_vertices, SubTreeReverse);
    
    int best_edge_disconnects_R = max_bridge_size;
    
    // swap adjacency lists
    temp_adj = First;
    First = r_First;
    r_First = temp_adj;
    
    temp_adj = Target;
    Target = r_Target;
    r_Target = temp_adj;
    
    temp_adj = Next;
    Next = r_Next;
    r_Next = temp_adj;
    
    temp_adj = Before;
    Before = r_Before;
    r_Before = temp_adj;
    
    if(best_edge_disconnects_F > best_edge_disconnects_R)	initReachableToInt(SubTreeReverse,0);
    else							initReachableToInt(SubTreeForward,0);

//     if(BridgesNum == 0){ //G is 2 edge connected
// 
//         if(temp_stack_top > max2vcbSize) max2vcbSize = temp_stack_top;max_bridge
//         if(temp_stack_top < min2vcbSize) min2vcbSize = temp_stack_top;
//         sum2vcbSize += temp_stack_top;
//         numOf2vcbs++; // number of 2ECC  
// 
//         //SortNumbers(temp_stack, temp_stack_top, n);
// 
//         printf("\n2ecc (%d vertices) : ",temp_stack_top);
//         for(int i = 0; i < temp_stack_top; i++) printf(" %d",temp_stack[i]);
//         printf("\n");
//     
//     }
    if(BridgesNum > 0){
        for(int i = 0; i < temp_stack_top; i++){
            int next_to_process = temp_stack[i];

            // remove the edges where the endpoint was marked in the previous step
            int id_forward = SubTreeForward[next_to_process];
            int id_reverse = SubTreeReverse[next_to_process];
	    
	    
	    int j = First[next_to_process];
	    while(j!=0){
		int v = Target[j];
		if(SubTreeForward[v] != id_forward || SubTreeReverse[v] != id_reverse){ 
		    delete_arc_from_graph(next_to_process,v,j);
		    //printf("[Bridges] edge (%d,%d) deleted\n",next_to_process,v);
		    stack_Cond[v]=true;
		    stack_Cond[next_to_process]=true;
		}
		j=Next[j];
	    }
        }
    }
    return BridgesNum;
}

//////// FIND 1EDGE OUT ////////




/* depth-first search from vertex k for 2*delta+1 edges 
   stores parents and labels        */

void DFS_1(int u, int *Target, int *Next, int *First, bool *marked_S1_Cond) {
    label2pre[u] = (++next);
    pre2label[next] = u;
   
    int k=First[u];
    while(k!=0){
	int v=Target[k];
	
	if (F_one >=(check)){
	    break;
	}
	
	F_one=F_one+1;
	w[u]=w[u]+1;
     
// 	printf("\n\n check (%d,%d)", u,v);
// 	printf("\n marked S1=%d",marked_S1[p]);
// 	printf("\n F_one=%d",F_one);    
// 	printf("\n check=%d",check);   
     
	if (!label2pre[v]) { 
	      marked_S1[p_S1]=v;
	      marked_S1_Cond[v]= true;
	      p_S1++;
	      DFS_1(v, Target, Next, First, marked_S1_Cond);
	      parent[v] = u;
	}
        k=Next[k];
    }  
}







/* depth-first search from vertex k for delta+1 edges 
   stores parents and labels        */

void DFS_2(int u, int *Target, int *Next, int *First, bool *marked_S2_Cond) {
    label2pre_2[u] = (++next_2);
   
    ///// first check for the reversed edge 
//     if (F_two > (delta)){
// 	printf("\n dfs has visited more than delta edges  \n ");
//      }
    
    if (w[u]>=delta && w[parent[u]]>=delta && F_two <= (delta)) {
	  F_two=F_two+1;
	  int v=parent[u];

	  if (!label2pre_2[v]) {
		marked_S2[p_S2]=v;
		marked_S2_Cond[v]= true;
		p_S2++;
		DFS_2(v, Target, Next, First, marked_S2_Cond);
	  }
    }
    
    int k = First[u];
    while(k!=0){
	int v=Target[k];
	if (F_two > (delta)){
	      break;
	}
	F_two=F_two+1;
	//printf("check DFS_2 (%d,%d)\n", u,v);
	/*printf("\n marked S1=%d",marked_S1[p]);
	printf("\n F_one=%d",F_two);    
	printf("\n check=%d",check); */  
//      
	if (!label2pre_2[v] && (u != parent[v] || w[v] < delta)) {
		marked_S2[p_S2]=v ;
		marked_S2_Cond[v]= true;
		p_S2++;
		DFS_2(v, Target, Next, First, marked_S2_Cond);  
	}
        k=Next[k];
    }  
} 	




int outgoing_edge_1;
int outgoing_edge_2;

int OneEdgeOut(int u, int *Target, int *Next, int *First, bool *marked_S1_Cond, bool *marked_S2_Cond){
    next=0;
    F_one=0;
    marked_S1[1]=u;
    marked_S1_Cond[u] = true;
    p_S1=2;
    parent[u]=0;
 	
    DFS_1(u, Target, Next, First, marked_S1_Cond);
 
    //printf("DFS 1 visited %d edges from %d | delta = %d | check = %d\n",F_one,u,delta,check);
    if (F_one < check){
      number_of_nodes_of_1egeOut=(p_S1);
      p=p_S1;
      idf=1;
      //printf("\n p=%d ,next=%d \n ",p,next);
      
      
      //printf("found component of type 1 from %d in G\n",u);
      //for(int i = 1; i < p_S1 ; i++) printf("Vertex %d belongs to the component\n",marked_S1[i]);
      return idf;
    }
    else{
		for (int i=next-1;i>1;i--){
		  //printf("lock (%d,%d)\n",parent[pre2label[i]],pre2label[i]);
		  w[parent[pre2label[i]]]+=w[pre2label[i]];
		}
		
		
		next_2=0;
		F_two=0;
		p_S2=2;
		marked_S2[1]=u;
		marked_S2_Cond[u] = true;
		
		DFS_2(u, Target, Next, First, marked_S2_Cond);
		
		//printf("DFS 2 visited %d edges from %d | delta = %d | check = %d\n",F_two,u,delta,check);
		if (F_two  < (delta+1)){
			//printf("found component of type 2 from %d in G\n",u);
			//for(int i = 1; i < p_S2 ; i++) printf("Vertex %d belongs to the component\n",marked_S2[i]);
			number_of_nodes_of_1egeOut=p_S2; 
			p=p_S2;
			idf=2;
			for(int i=next-1;i>1;i--){
				if(w[pre2label[i]]>delta+1 && !marked_S2_Cond[pre2label[i]] && marked_S2_Cond[parent[pre2label[i]]]){
					outgoing_edge_1 = parent[pre2label[i]];
					outgoing_edge_2 = pre2label[i];
					break;
				}
			}
			return idf;
		}
		else{   
			number_of_nodes_of_1egeOut=0;
			return 0;
		}
    }
}

  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////check for one-edges-out-component ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Process_OneEdgeOut(int vertex, int size_of_SCC, int num_of_edges, int recursionDepth){
	
	if(recursionDepth > maxRecursionDepth) maxRecursionDepth = recursionDepth;
	numOfRecursionCalls ++;
	
	//printf("Process_OneEdgeOut(%d,%d,%d, %d)\n",vertex,size_of_SCC, num_of_edges,recursionDepth);
	int component_type;
	int component_type_r;
	
	int v,l,k;
	//printf("\n stack top= %d \n",stack_top);
	
	fillReachableStack(vertex, temp_stack);
	
	int num_of_SBs = deleteBridges(vertex, size_of_SCC);
	
	if(num_of_SBs == 0){
		if(temp_stack_top > max2ECSSize) max2ECSSize = temp_stack_top;
		if(temp_stack_top < min2ECSSize) min2ECSSize = temp_stack_top;
		sum2ECSSize += temp_stack_top;
		numOf2ECS++; // number of 2ECC  
		
		//SortNumbers(temp_stack, temp_stack_top, n);
		
		//         printf("\n2ecc (%d vertices) : ",temp_stack_top);
		//         for(int i = 0; i < temp_stack_top; i++)// printf(" %d",temp_stack[i]);
		//         printf("\n");
		return;
	}
	
	for (int i=0;i<size_of_SCC+2;i++){
		pre2label[i]=0;
	}
	
	initReachableToInt(label2pre, 0);
	initReachableToInt(label2pre_2, 0);
	initReachableToInt(marked_S1, 0);
	initReachableToInt(marked_S2, 0);
	initReachableToInt(w, 0);
	initReachableToFalse(marked_Cond);
	initReachableToFalse(marked_S1_Cond);
	initReachableToFalse(marked_S2_Cond);
	/*
	for(int i=0 ; i <= n; i++){
	label2pre[i] = 0;
	label2pre_2[i] = 0;
	marked_S1[i] = 0;
	marked_S2[i] = 0;
	w[i] = 0;
	marked_Cond[i] = false;
	marked_S1_Cond[i] = false;
	marked_S2_Cond[i] = false;
	}*/
	
	stack_top = 0;
	for (int i = 0; i < temp_stack_top; i++) {
	    if (stack_Cond[temp_stack[i]]){
		stack[++stack_top] = temp_stack[i]; 
	    }
	}

	while (stack_top > 0) {
		int u = stack[stack_top--];
		
		if (stack_Cond[u] == false || in_small_comp[u]) continue;
		//printf("\nstack top (%d)= %d \n",stack_top,u);
		stack_Cond[u] = false;
		//int u = stack[stack_top--];
		//printf("stack top (%d)= %d \n",prev_stack_of_u,u);
		number_of_nodes_of_1egeOut=0;
		p_S1 = p_S2 = 0;
		
		component_type = OneEdgeOut(u, Target, Next, First, marked_S1_Cond, marked_S2_Cond);
		
		//printf("number_of_nodes_of_1egeOut= %d, check=%d\n ",number_of_nodes_of_1egeOut,check);
		
		if (component_type != 0){
			number_of_components = number_of_components+1;
			
			
			if(number_of_nodes_of_1egeOut > max1EOutSize) max1EOutSize = number_of_nodes_of_1egeOut;
			if(number_of_nodes_of_1egeOut < min1EOutSize) min1EOutSize = number_of_nodes_of_1egeOut;
			sum1EOutSize += number_of_nodes_of_1egeOut;
			numOf1EOut++; // number of 2ECC  
			//if(number_of_nodes_of_1egeOut> 40) printf("\n number_of_nodes_of_1egeOut= %d, check=%d \n ",number_of_nodes_of_1egeOut,check);
			
			
			int *S_one;
			if(idf == 2){
				l = First[outgoing_edge_1];
				while(l!=0){
					if(Target[l]==outgoing_edge_2){
						delete_arc_from_graph(outgoing_edge_1,outgoing_edge_2,l);
						//printf("1)-----> edge (%d,%d) deleted\n",outgoing_edge_1,outgoing_edge_2);
						
						if (!in_small_comp[outgoing_edge_2] && stack_Cond[outgoing_edge_2] == false) {
							stack[++stack_top] = outgoing_edge_2;
							stack_Cond[outgoing_edge_2] = true;
						}
					}
					l=Next[l];
				}
				S_one= marked_S2;
			}
			else{
				S_one= marked_S1;
			}
			
			/// delete inicident edges from the graph 
			for (int i=1;i<number_of_nodes_of_1egeOut;i++){
			  
				l=r_First[S_one[i]];
				while(l!=0){
					v=r_Target[l]; 
					if (idf==1){
						if(marked_S1_Cond[v]==true){
							//printf("\n the edge belongs to the component\n");
							l=r_Next[l];
							continue;
						}
					}
					else if (idf==2){
						if (marked_S2_Cond[v]==true){
							// printf("\n the edge belongs to the component\n");
							l=r_Next[l];
							continue;
						}
					}
					
					//printf("\n  delete incident edge (%d,%d) from forward graph \n",v, S_one[i]);
					//printf("|");	
					delete_arc_from_graph(v,S_one[i],l);
					//printf("2)-----> edge (%d,%d) deleted\n",v,S_one[i]);
					
					if (!in_small_comp[v] && stack_Cond[v] == false) {
						stack[++stack_top] = v;
						stack_Cond[v] = true;
					}
					stack_Cond[S_one[i]] = false;
					in_small_comp[S_one[i]] = true;
					
					l = r_Next[l];
				}
			}
		}
		
		//// initialize the structures that used in OneEdgeOut
		if (p_S1>p_S2){
			k=p_S1;
		}
		else{
			k=p_S2;
		}
		
		
		for (int i=1;i<=k;i++){
			marked_S1_Cond[marked_S1[i]]=false;
			marked_S2_Cond[marked_S2[i]]=false;
			
			pre2label[i]=0;
			
			label2pre[marked_S1[i]]=0;
			label2pre_2[marked_S2[i]]=0;
			
			w[marked_S1[i]]=0;
			
			marked_S1[i]=0;
			marked_S2[i]=0;
		}
		
		if(component_type == 0){
			//////////////////////////////////////////////  
			number_of_nodes_of_1egeOut=0;
			//////////////  REVERSE GRAPH //////////////////////////
			
			p_S1 = p_S2 = 0;
			component_type_r = OneEdgeOut(u, r_Target, r_Next, r_First, marked_S1_Cond, marked_S2_Cond);
			
			if (component_type_r != 0){
				//printf("found component from %d in G^R\n",u);
			  
				if(number_of_nodes_of_1egeOut > max1EOutSize) max1EOutSize = number_of_nodes_of_1egeOut;
				if(number_of_nodes_of_1egeOut < min1EOutSize) min1EOutSize = number_of_nodes_of_1egeOut;
				sum1EOutSize += number_of_nodes_of_1egeOut;
				numOf1EOut++; // number of 2ECC  
				number_of_components=number_of_components+1;
				//if(number_of_nodes_of_1egeOut> 40) printf("\n number_of_nodes_of_1egeOut= %d, check=%d \n ",number_of_nodes_of_1egeOut,check);
				/// delete inicident edges from the graph 
				
				int *S_one_r;
				if(idf == 2){
					l = r_First[outgoing_edge_1];
					while(l!=0){
						if(r_Target[l]==outgoing_edge_2){
							delete_arc_from_graph(outgoing_edge_2,outgoing_edge_1,l); 
							//printf("3)-----> edge (%d,%d) deleted\n",outgoing_edge_2,outgoing_edge_1);
							if (!in_small_comp[outgoing_edge_2] && stack_Cond[outgoing_edge_2] == false) {
								stack[++stack_top] = outgoing_edge_2;
								stack_Cond[outgoing_edge_2] = true;
							}
						}
						l=r_Next[l];
					}
					S_one_r = marked_S2;
				}
				else{
					S_one_r = marked_S1;
				}
				
				//else{continue;}
				for (int i=1;i<number_of_nodes_of_1egeOut;i++){
					l=First[S_one_r[i]];
					while(l!=0){
						v=Target[l]; 
						
						if (idf==1){
							if (marked_S1_Cond[v]==true){
								// printf("\n the edge belongs to the component\n");
								l=Next[l];
								continue;
							}
						}
						else if (idf==2){
							if (marked_S2_Cond[v]==true){
								// printf("\n the edge belongs to the component\n");
								l=Next[l];
								continue;
							}
						}  
						
						//printf("\n  delete incident edge (%d,%d) from reverse graph \n",v, S_one[i]);
						
						delete_arc_from_graph(S_one_r[i],v,l);	

						if (!in_small_comp[v]) {
						    //add their endpoints in stack
						    if (stack_Cond[v] == false) {
							stack[++stack_top] = v;
							stack_Cond[v] = true;
						    }
						}

						stack_Cond[S_one_r[i]] = false;
						in_small_comp[S_one_r[i]] = true;
			    
						l = Next[l];
					}
				}
			}
			
			//// initialize the structures that used in OneEdgeOut
			if (p_S1>p_S2){
				k=p_S1;
			}
			else{
				k=p_S2;
			}
			
			
			for (int i=1;i<=k;i++){
				marked_S1_Cond[marked_S1[i]]=false;
				marked_S2_Cond[marked_S2[i]]=false;
				
				pre2label[i]=0;
				
				label2pre[marked_S1[i]]=0;
				label2pre_2[marked_S2[i]]=0;
				
				w[marked_S1[i]]=0;
				
				marked_S1[i]=0;
				marked_S2[i]=0;
			}
		}
	}
	// compute the SCCs of the remaining graph
	int currentSccID = 1;
	connectedComponents(&currentSccID);
	
	int *num_in_SCC = new int [currentSccID];
	int *vertex_in_SCC = new int [currentSccID];
	
	int x;
	for(int i = 1; i < currentSccID; i++){
		for(int j = SCCs_first[i]; j < SCCs_first[i+1]; j++){ // for every vertex in the SCC
			x = SCCs[j];
			l=First[x];
			while(l!=0){
				v=Target[l]; 
				if(SCC_ID[v] != i){
					delete_arc_from_graph(x,v,l);
					stack_Cond[v]=true;
					stack_Cond[x]=true;
				}
				l=Next[l];
			}
		}
	}
	
	
	int *swap_reachable  = new int [temp_stack_top];
	int swap_reachable_top = temp_stack_top;
	for(int i = 0; i < temp_stack_top; i++) swap_reachable[i] = temp_stack[i];
	
	// store the size and a vertex for each SCC (because the SCC matrices are global)
	for(int i = 1; i < currentSccID; i++){
		num_in_SCC[i] = SCCs_first[i+1] - SCCs_first[i];
		vertex_in_SCC[i] = SCCs[SCCs_first[i]];
	}
	
	for(int i = 1; i < currentSccID; i++){
		if(num_in_SCC[i] > 2){
			// recursive call
			int num_of_SBs = deleteBridges(vertex_in_SCC[i], num_in_SCC[i]);
		}
	}
	
	for(int i = 0; i < swap_reachable_top; i++) temp_stack[i] = swap_reachable[i];
	temp_stack_top = swap_reachable_top;
	delete [] swap_reachable;
	delete [] num_in_SCC;
	delete [] vertex_in_SCC;
	
	
// 	printf("%d SCCs were found\n",currentSccID-1);
// 	exit(1);
	
	currentSccID = 1;
	connectedComponents(&currentSccID);
	
	num_in_SCC = new int [currentSccID];
	int *edges_in_SCC = new int [currentSccID];
	vertex_in_SCC = new int [currentSccID];
	
	for(int i = 1; i < currentSccID; i++){
		edges_in_SCC[i] = 0;
	}
	
	for(int i = 1; i < currentSccID; i++){
		for(int j = SCCs_first[i]; j < SCCs_first[i+1]; j++){ // for every vertex in the SCC
			x = SCCs[j];
			l=First[x];
			while(l!=0){
				v=Target[l]; 
				if(SCC_ID[v] != i){
					delete_arc_from_graph(x,v,l);
					stack_Cond[v]=true;
					stack_Cond[x]=true;
				}
				else{
					edges_in_SCC[i]++;
				}
				l=Next[l];
			}
		}
	}
	
	// store the size and a vertex for each SCC (because the SCC matrices are global)
	for(int i = 1; i < currentSccID; i++){
		num_in_SCC[i] = SCCs_first[i+1] - SCCs_first[i];
		vertex_in_SCC[i] = SCCs[SCCs_first[i]];
	}
	
	for(int i = 1; i < currentSccID; i++){
		if(num_in_SCC[i] > 2){
			//recursive call
			// 	    printf("Call %d-th SCC (out of %d)\n",i,currentSccID-1);
			Process_OneEdgeOut(vertex_in_SCC[i], num_in_SCC[i], edges_in_SCC[i],recursionDepth+1);
		}
	}
	
	delete [] num_in_SCC;
	delete [] vertex_in_SCC;
	delete [] edges_in_SCC;
}

 

void delete_space () {
     delete [] First;
     delete [] Target;
     delete [] Next;
     delete [] Before;
     delete [] r_First;
     delete [] r_Target;
     delete [] r_Next;
     delete [] r_Before;
     delete [] label2pre; 
     delete [] pre2label ;
     delete [] parent ;
     delete []label2pre_2; 
     delete [] marked_Cond ;
   
     delete [] marked_S1_Cond;
     delete [] marked_S2_Cond;
    
     delete [] stack_Cond;
     delete [] in_small_comp;
    
     delete [] w;
    
     delete [] LowLink ;
    
     delete [] marked_S1;
     delete [] marked_S2;
     delete [] stack ;
     delete [] temp_stack ;
     delete [] temp_inserted ;
    
     delete [] SCCs;
     delete [] SCCs_first;
     delete [] SCC_ID ;
  
  
}
void init() {
    current_pos = r_current_pos = 1;
  
    First=new int [n+2];
    Target=new int [m+2];
    Next=new int [m+2];
    Before=new int [m+2];
    
    
    r_First=new int [n+2];
    r_Target=new int [m+2];
    r_Next=new int [m+2];
    r_Before=new int [m+2];
    label2pre = new int [n+2];
    pre2label = new int [n+2];
    parent = new int [n+2];
    label2pre_2 = new int [n+2];
    marked_Cond = new bool [n+2];
   
    marked_S1_Cond=new bool [n+2];
    marked_S2_Cond=new bool [n+2];
    
    stack_Cond=new bool [n+2];
    in_small_comp=new bool [n+2];
    
    w = new int [n+2];
    
    LowLink = new int [n+1];
    
    marked_S1=new int [n+2];
    marked_S2=new int [n+2];
    stack = new int [2*n+2];
    temp_stack = new int [n+2];
    temp_inserted = new bool[n+2];
    
    SCCs = new int [n+2];
    SCCs_first= new int [n+2];
    SCC_ID = new int [n+2];
    
    max2ECSSize = sum2ECSSize = numOf2ECS = maxRecursionDepth = numOfRecursionCalls = 0;
    max1EOutSize=sum1EOutSize=numOf1EOut =0;
    min2ECSSize=min1EOutSize=n;

    for (int i=0;i<n+2;i++){
      
		stack[i]=i;
		
	// 	label2pre[i]=0;
	// 	pre2label[i]=0;
	// 	
	// 	label2pre_2[i]=0;
	// 	pre2label_2[i]=0;
	// 	r_label2pre[i]=0; 
	// 	r_pre2label[i]=0;
		
		First[i]=0; 
		r_First[i]=0;
		
		
		in_small_comp[i] = false;
		stack_Cond[i] = true; // here it might as well be set to true, which would mean that in the beginning of the algorithm there will be searches from each vertex
		temp_inserted[i] = false;
    }
   
    for (int i=0;i<m+2;i++){
      Next[i]=0;
      Before[i]=0;
      r_Next[i]=0;
      r_Before[i]=0;
    }
}


void builtgraph() {
    int input_source,input_target;
    
    for (int current_input_pos=0 ; current_input_pos < m; current_input_pos++){
        input_source=input_file[2*current_input_pos];
        input_target=input_file[2*current_input_pos+1];
		//printf("\n insert (%d,%d)", input_source,input_target );
		insert_arc_to_graph(input_source,input_target);
	
     }
}



/* read graph from input file */
void readGraph(const char *file) {
    FILE *input = fopen(file, "r");
    if (!input) {
        fprintf(stderr, "Error opening file \"%s\".\n", file);
        exit(-1);
    }

    int x, y;
    int p = 0;

    while (fgets(line, MAXLINE, input) != NULL) {
        switch (line[0]) {
            case '\n':;
            case '\0': break;
            case 'p': if (sscanf(line, "p %d %d", &n, &m) != 2) {
                    fprintf(stderr, "Error reading graph size (%s).\n", file);
                    exit(-1);
                }
                input_file = new int [2*m];
                break;
            case 'a': if (sscanf(line, "a %d %d", &x, &y) != 2) {
                    fprintf(stderr, "Error reading graph arc (%s).\n", file);
                    exit(-1);
                }
                assert(x <= n);
                assert(y <= n);

                if(x == y) break;

                input_file[p++] = x;
                input_file[p++] = y;
                if (p>2*m) {
                    fprintf(stderr, "Error! Graph has >%d arcs.\n", m);
                    exit(-1);
                }
                break;
            default: fprintf(stderr, "Error reading graph (%s).\n", file);
                exit(-1);
        }
    }
    printf("number of nodes =%d , number of edges=%d\n",n,m);
    fprintf(stderr, "END reading graph (%s).\n", file);
    fclose(input);
}




 int main(int argc, char *argv[]) {
    if (argc != 4 ) {
        printf("\n usage: %s <input file> <seed> <starting vertex> \n\n", argv[0]);
        exit(-1);
    }

    char *file = argv[1];
    int seed = atoi(argv[2]);
    int rand_node = atoi(argv[3]);
    //delta= atoi(argv[4]);
    //int prm= atoi(argv[4]);

    readGraph(file);
  
    
    //To Create the output File
    inputFileName = Utilities::TO_STRING(file);
    preText = Utilities::GetFileName(inputFileName);
    testDate = Utilities::GetCurrentDateTime();

    RFWTimer timer(true);
    
    init();
    builtgraph();
 
    
    srand(seed);
    delta = sqrt(m);
    check = 2 * delta + 1; // check=2*delta +1
    printf("\n delta =%d \n ",delta);
    Process_OneEdgeOut(rand_node,n,m,1);
    delete_space();

    
    takenTime = timer.getTime();

    if (numOf2ECS > 0) {
        avg2eccSize = (float) sum2ECSSize / (float) numOf2ECS;
    }

     ShowOutput();
     if (is_export_output)
         ExportOutput();

     // (For the experiment, we don't need the following codes)
     return 0;
    
    if(numOf2ECS > 0){
        fprintf (stdout, "Number of 2ECCs       : %d\n", numOf2ECS);
        fprintf (stdout, "Average 2ECC size     : %f\n", (float)sum2ECSSize / (float)numOf2ECS);
        fprintf (stdout, "Max 2ECC              : %d\n", max2ECSSize);
        fprintf (stdout, "Min 2ECC              : %d\n", min2ECSSize);
    }
    
    fprintf (stdout,     "====================================\n");
    fprintf (stdout,     "Total calls to main		: %d\n", numOfRecursionCalls);
    fprintf (stdout,     "Max recursion depth		: %d\n", maxRecursionDepth);
    
    fprintf (stdout,     "====================================\n");
    //fprintf (stdout, "Delta = %d, 2*Delta+1 = %d       : %d\n", delta, check);
    fprintf (stdout, "Number of 1EOut components       : %d\n", numOf1EOut);
    fprintf (stdout, "Average 1EOut component size     : %f\n", (float)sum1EOutSize / (float)numOf1EOut);
    fprintf (stdout, "Max 1EOut component              : %d\n", max1EOutSize);
    fprintf (stdout, "Min 1EOut component              : %d\n", min1EOutSize);
    
    cout << " Total Time Taken         : " << takenTime << " seconds"
            << endl;
    //fprintf (stdout, "totaltime %f\n", t);
   
}


void ShowOutput() {

    // Source code added by Nilu

    cout << " x==========================================x" << endl;
    cout << " |     Summary of the 2ESCC Computation     |" << endl;
    cout << " x==========================================x" << endl;

    // Algorithm code, test date and input file
    cout << " Algorithm-Code           : " << algorithm_code << endl;
    cout << " Test Date                : " << testDate << endl;
    cout << " Input Graph File         : " << inputFileName << endl;

    // Graph Details and Ratio
    cout << " Graph Details G(V,E)     : [" << n << "," << m << "]" << endl;

    cout << " Ratio (E:V)              : " << Utilities::GetRatioString(n, m)
            << endl;

    // Computational Result
    cout << " Total 2ESCC components   : " << numOf2ECS << endl;
    cout << " Total Time Taken         : " << takenTime << " seconds" << endl;
    cout << " Recursion Depth Level    : " << maxRecursionDepth << endl;
    cout << " Minimum 2ESCC Vertex Size: " << min2ECSSize << endl;
    cout << " Maximum 2ESCC Vertex Size: " << max2ECSSize << endl;
    cout << " Average 2ESCC Vertex Size: " << avg2eccSize << endl;
    cout << " Total Fractions          : " << numOfRecursionCalls << endl;

}

void ExportOutput() {

    string fileName = preText + "_2ESCC_TIME_" + algorithm_code + ".txt";
    ofstream out_file;
    out_file.open(fileName.c_str());

    //Header on file
    out_file << " x==========================================x" << endl;
    out_file << " |     Summary of the 2ESCC Computation     |" << endl;
    out_file << " x==========================================x" << endl;

    // Algorithm code, test date and input file
    out_file << " Algorithm-Code           : " << algorithm_code << endl;
    out_file << " Test Date                : " << testDate << endl;
    out_file << " Input Graph File         : " << inputFileName << endl;

    // Graph Details and Ratio
    out_file << " Graph Details G(V,E)     : [" << n << "," << m << "]" << endl;

    out_file << " Ratio (E:V)              : "
            << Utilities::GetRatioString(n, m) << endl;

    // Computational Result
    out_file << " Total 2ESCC components   : " << numOf2ECS << endl;
    out_file << " Total Time Taken         : " << takenTime << " seconds"
            << endl;
    out_file << " Total Recursion Depth    : " << maxRecursionDepth << endl;

    out_file << " Minimum 2ESCC Vertex Size: " << min2ECSSize << endl;
    out_file << " Maximum 2ESCC Vertex Size: " << max2ECSSize << endl;
    out_file << " Average 2ESCC Vertex Size: " << avg2eccSize << endl;
    out_file << " Total Fractions          : " << numOfRecursionCalls << endl;

    out_file.close();

}
