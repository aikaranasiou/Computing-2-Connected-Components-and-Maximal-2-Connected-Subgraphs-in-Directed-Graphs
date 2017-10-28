#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <string.h>
#include "rfw_timer.h"
#include <math.h>

#include "utilities.h"

// using namespace std;

#define MAXLINE       100   /* max length of input line */
char line[MAXLINE]; /* stores input line */

int Nodes, M;  /* number of nodes, arcs */
int n, m; /* double number of nodes ,arcs */ 
int *input_file; // stores input arcs
FILE *output;
int x_sp;
int j_global=0;

// dynamic list
int *First;
int *Target;
int *Next;
int *Before;

int *r_First;
int *r_Target;
int *r_Next;
int *r_Before;

int *real_arc;
int *marked_S1;
int *marked_S2;
int starting_vertex;
int *LowLink;
bool *SAP_in_Forward;
bool *SAP_in_Reverse;

int *Nearest_SAP_Forward;
int *Nearest_SAP_Reverse;
int SAP;

int *stack ;

int *marked_BFS ;
int *SCCs ;
int *SCCs_first ;
int *SCC_ID ;

int *queue_BFS;
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
int stack_top; 
bool *stack_Cond;
bool *marked_S1_Cond;
bool *marked_S2_Cond;
int *w_matrix;

int maxRecursionDepth, max2VCSSize, min2VCSSize, sum2VCSSize, numOf2VCS, numOfRecursionCalls;
int max1VOutSize, min1VOutSize, sum1VOutSize, numOf1VOut;


int *temp_stack, temp_stack_bottom, temp_stack_top;
bool *temp_inserted;
int k_split;

int next_free;

bool is_export_output = true;
string algorithm_code = "2V-CHILP";
float avg2eccSize = 0;
string inputFileName = "";
string preText = "";
string testDate = "";
double takenTime = 0;

void ShowOutput();
void ExportOutput();

inline void insert_arc_to_graph(int x,int y){
    int future_next_free = Next[next_free];
    if(First[x]==0){
        Target[next_free]=y;
        First[x]=next_free;
        Before[next_free] = 0;
	Next[next_free] = 0;
    }
    else{
        Target[next_free]=y;
        Next[next_free]=First[x];
	Before[First[x]]=next_free;
        First[x]=next_free;
    }
    if(r_First[y]==0){
        r_Target[next_free]=x;
        r_First[y]=next_free;
        r_Before[next_free] = 0;
	r_Next[next_free] = 0;
    }
    else{
        r_Target[next_free]=x;
        r_Next[next_free]=r_First[y];
	r_Before[r_First[y]]=next_free;
        r_First[y]=next_free;
    }
    next_free = future_next_free;
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
    Next[pos] = next_free;
    next_free = pos;
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
	      int future_next = Next[j];
	      if(!temp_inserted[v]){
		  stack[temp_stack_top++] = v;
		  temp_inserted[v] = true;
	      }
	      j=future_next;
	}
	
	j = r_First[temp_vertex];
	while(j!=0){
	      int v = r_Target[j];
	      int future_next = r_Next[j];
	      if(!temp_inserted[v]){
		  stack[temp_stack_top++] = v;
		  temp_inserted[v] = true;
	      }
	      j=future_next;
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



inline void rcompress(int v, int *parent, int *semi, int *label, int c) {
    int p;
    if ((p = parent[v]) > c) {
        rcompress(p, parent, semi, label, c);
        if (semi[label[p]] < semi[label[v]]) label[v] = label[p];
        parent[v] = parent[p];
    }
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
	int future_next = Next[j];
	if ( !label2pre[v] ){
            dfs_parent[v] = -k; // use negative for vertices that have not been processed
            DFS(v, N, parent, skip);
            parent[label2pre[v]] = label2pre[k];
        }
        else if ( dfs_parent[v] > 0 ) skip[v] = true; // (k,v) is forward or cross arc
        
	j=future_next;
    }

    dfs_parent[k] = -dfs_parent[k]; // k is processed
}

void StrongConnect(int k, int *N, int *currentSccID, int *SCC_next_pos, bool* isInStack){

    int q;
    LowLink[k] = label2pre[k] = ++(*N);
    stack[++stack_top] = k;
    isInStack[k] = true;

    int j = First[k];
    while(j!=0){
	  int v = Target[j];
	  int future_next = Next[j];
	  
	  if(!label2pre[v]){ // tree arc
	      StrongConnect(v, N, currentSccID, SCC_next_pos, isInStack);
	      LowLink[k] = min(LowLink[k], LowLink[v]);
	      
	  }
	  else if(label2pre[v] < label2pre[k]){ // front or cross arc
	      if(isInStack[v]) LowLink[k] = min(LowLink[k], label2pre[v]);
	  }
	  j=future_next;
    }
    if(LowLink[k] == label2pre[k]){ // k is the root of a component
        while( (stack_top > 0) && (label2pre[q = stack[stack_top]] >= label2pre[k]) ){
            SCC_ID[q] = (*currentSccID);
            isInStack[q] = false;
            SCCs[(*SCC_next_pos)++] = q;
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
    stack_top = 0;

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




/* depth-first search of dominator tree from node k; assigns preorder numbers and sizes */
void DomDFS(int k, int *dN, int *d_adj, int *d_adj_first, int *dlabel2pre, int *dsize)
{
    dlabel2pre[pre2label[k]]  = ++(*dN);
    dsize[pre2label[k]] = 1;

    for (int j=d_adj_first[k]; j<d_adj_first[k+1]; j++) {
        int v = d_adj[j];
        DomDFS(v, dN, d_adj, d_adj_first, dlabel2pre, dsize);
        dsize[pre2label[k]] += dsize[pre2label[v]];
    }
    //printf("dsize[%d]=%d\n",pre2label[k],dsize[pre2label[k]]);
}

void DomDFS_mark(int k, int *d_adj, int *d_adj_first, bool *isSap, int *Nearest_SAP)
{
    if(isSap[pre2label[k]]) Nearest_SAP[pre2label[k]] = pre2label[k];
    
    for (int j=d_adj_first[k]; j<d_adj_first[k+1]; j++) {
        int v = d_adj[j];
	Nearest_SAP[pre2label[v]] = pre2label[k];
        DomDFS_mark(v, d_adj, d_adj_first, isSap, Nearest_SAP);
    }
}

/* constructs the dominator tree as a graph and perfors a DFS*/
inline void ConstructDomTree(int verticesNum, int *Dout, int *DfirstOut, int *dom) {
    int bsize = verticesNum + 2;
    int *DnextOut  = new int [bsize];

    for (int i = 0; i <= verticesNum; i++) {
        DfirstOut[i] = 0;
    }
 
    DfirstOut[verticesNum+1] = DnextOut[verticesNum+1] = 0;
    
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

/*computed dominators of a flow graph rooted at r1*/
int slt(int r1, int num_of_vertices, int *idom, bool *isSap, int *Nearest_SAP, int *dlabel2pre, int *dsize) {

    int bsize = num_of_vertices + 1;
    int *parent = new int [bsize];
    int *semi = new int [bsize];
    int *label= new int [bsize];
    int *ubucket = new int [bsize];
    int *dom = new int [bsize];
    int i;

    dfs_parent = marked_S1;
    
    for (i = num_of_vertices; i >= 0; i--) {
        label[i] = semi[i] = i;
        ubucket[i] = 0;
        dom[i] = 0;
    }

    // initialize matrices only for the reachable vertices
    //initReachableToInt(dfs_parent, 0);
    fillReachableStack(r1, temp_stack);
    initReachableToInt(label2pre, 0);
    
    
    bool *skip = marked_Cond;
    initReachableToFalse(skip);
    initReachableToFalse(isSap);

    int N = 0;

    parent[1] = 1;
    DFS(r1, &N, parent, skip);
    
//      printf("visited %d/%d vertices\n",N,num_of_vertices);

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
	    int future_next = r_Next[j];
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
	    j=future_next;
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
    idom[pre2label[1]] = pre2label[dom[1]];
    for ( i = 2; i <= N; i++) {
        if (dom[i] != semi[i]) dom[i] = dom[dom[i]]; //make relative absolute // stores preorder
        idom[pre2label[i]] = pre2label[dom[i]];
// 	printf("idom[%d]=%d\n",pre2label[i],idom[pre2label[i]]);
    }
    
    
    int *Dout      = new int [N+2];
    int *DfirstOut = new int [N+2];
    
    // construct dominator tree
    ConstructDomTree(N, Dout, DfirstOut, dom);
    
    int dN = 0;
    DomDFS(1, &dN, Dout, DfirstOut, dlabel2pre, dsize);
    
   int num_of_SAPs = 0;
   int best_split_vertex = 0;
   int best_split = 0;
    for (i = 2; i <= N; i++) {
      if (dom[i]!=1){
	  if(dsize[pre2label[dom[i]]] > best_split && (num_of_vertices-dsize[pre2label[dom[i]]])> best_split){
	    if(dsize[pre2label[dom[i]]] > (num_of_vertices-dsize[pre2label[dom[i]]])) best_split = (num_of_vertices-dsize[pre2label[dom[i]]]);
	    else best_split = dsize[pre2label[dom[i]]];
	    best_split_vertex = pre2label[dom[i]];
	  }
	  
// 	  isSap[pre2label[dom[i]]] = true;
// 	  num_of_SAPs = 1;
      }
    }
    if(best_split != 0){
      num_of_SAPs = 1;
      SAP = best_split_vertex;
      isSap[SAP] = true;
    }

    DomDFS_mark(1, Dout, DfirstOut, isSap, Nearest_SAP);
    
    delete [] Dout;
    delete [] DfirstOut;
    
    delete [] parent;
    delete [] semi;
    delete [] label;
    delete [] ubucket;
    delete [] dom;
    
    return num_of_SAPs;
}

/* depth-first search from vertex k excluding vertex x
   only keeps track of visited vertices                */
void xDFS(int k, int *N, int *First, int *Target, int *Next, int *label2pre, int x) {
    int j;

    label2pre[k] = ++(*N);

    j = First[k];
    while(j!=0){
	int v = Target[j];
	int future_next = Next[j];
	if ( !label2pre[v] && (v!=x)){
            xDFS(v, N, First, Target, Next, label2pre, x);
        }
        j=future_next;
    }
}

// test strong connectivity of G-r1
void testsc(int r1, int r2, int num_of_vertices){

    bool reach, reached;
    int N;

    reach = reached = false;
    // test if all vertices in V-r1 are reachable from r2
    initReachableToInt(label2pre,0);

    N = 0;
    xDFS(r2, &N, First, Target, Next, label2pre, r1);

    if (N != (num_of_vertices-1)) reach = false;
    else reach = true;

    // test if all vertices in V-r1 reach r2
    initReachableToInt(label2pre_2,0);
    N = 0;
    xDFS(r2, &N, r_First, r_Target, r_Next, label2pre_2, r1);
    if (N != (num_of_vertices-1)) reached = false;
    else reached = true;

    // return if the r1 is a strong articulation point
    if(!reach || !reached){
	j_global = j_global + 1;
	int x = Nodes + j_global;
	real_arc[x] = real_arc[r1];
	stack_Cond[x] = true;
	stack_Cond[r1] = true;
	
	int j=First[r1];
	while(j!=0){
	    int v = Target[j];
	    int future_next  = Next[j];
	    if (label2pre[v] && label2pre_2[v]){
	      //split
	      delete_arc_from_graph(r1,v,j);
	      //printf("delete edge (%d,%d) to descendant in D\n",SAP,v);
	      insert_arc_to_graph(x,v);
	    }
	    j=future_next;
	}
	
	j=r_First[r1];
	while(j!=0){
	    int v = r_Target[j];
	    int future_next  = r_Next[j];
	    if (label2pre[v] && label2pre_2[v]){
	      delete_arc_from_graph(v,r1,j);
	      //printf("delete edge (%d,%d) to descendant in D^R\n",SAP,v);
	      insert_arc_to_graph(v,x);
	    }
	    j=future_next;
	}
    }
}

int find_SAPS(int vertex, int size_of_SCC){
    int * dlabel2pre = label2pre_2;
    int * dsize = marked_S2;
    int * idom = w_matrix;
    Nearest_SAP_Forward = LowLink;
    int num_SAP_Forward = slt(vertex, size_of_SCC, idom, SAP_in_Forward, Nearest_SAP_Forward, dlabel2pre, dsize);
    
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
  
    int * r_dlabel2pre = w;
    int * r_dsize = stack;
    int * r_idom = marked_BFS;
    Nearest_SAP_Reverse = queue_BFS;
    int num_SAP_Reverse = slt(vertex, size_of_SCC, r_idom, SAP_in_Reverse, Nearest_SAP_Reverse, r_dlabel2pre, r_dsize);
    
    int num_SAP = num_SAP_Forward + num_SAP_Reverse;
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
    
	  
    int Nodes_before = Nodes+j_global;
    if(num_SAP > 0){
	//printf("vertex %d is a strong articulation point, min size in dom = %d\n",SAP, best_split);
	if(num_SAP_Forward > num_SAP_Reverse){
	  for(int i=0; i <temp_stack_top; i++){
	    int SAP = temp_stack[i];
	    if(SAP_in_Forward[SAP]){
		  j_global = j_global + 1;
		  int x = Nodes + j_global;
		  real_arc[x] = real_arc[SAP];
		  stack_Cond[x] = true;
		  stack_Cond[SAP] = true;
		  
		  int j=First[SAP];
		  while(j!=0){
		      int v = Target[j];
		      int future_next  = Next[j];
		      if(v > Nodes_before){
			j = future_next;
			continue;
		      }
		      if (idom[v]==SAP){
			//split
			delete_arc_from_graph(SAP,v,j);
			//printf("delete edge (%d,%d) to descendant in D\n",SAP,v);
			insert_arc_to_graph(x,v);
		      }
		      j=future_next;
		  }
		  
		  j=r_First[SAP];
		  while(j!=0){
		      int v = r_Target[j];
		      int future_next  = r_Next[j];
		      if(v > Nodes_before){
			j = future_next;
			continue;
		      }
		      if ( (dlabel2pre[v] >= dlabel2pre[SAP] ) && ( dlabel2pre[v] < dlabel2pre[SAP] + dsize[SAP] ) ) {
			//split
			delete_arc_from_graph(v,SAP,j);
			if(idom[v] == SAP || Nearest_SAP_Forward[v] == SAP){
			  insert_arc_to_graph(v,x);
			}
			else{
			  stack_Cond[v] = true;
			}
		      }
		      j=future_next;
		  }
	    }
	  }
	}
	else{
	  for(int i=0; i <temp_stack_top; i++){
	    int SAP = temp_stack[i];
	    if(SAP_in_Reverse[SAP]){
		j_global = j_global + 1;
		int x = Nodes + j_global;
		real_arc[x] = real_arc[SAP];
		stack_Cond[x] = true;
		stack_Cond[SAP] = true;
		
		int j = r_First[SAP];
		while(j!=0){
		    int v = r_Target[j];
		    int future_next  = r_Next[j];
		    if(v > Nodes_before){
		      j = future_next;
		      continue;
		    }
		    if (r_idom[v]==SAP){
		      //split
		      delete_arc_from_graph(v,SAP,j);
		      insert_arc_to_graph(v,x);
		    }
		    j=future_next;
		}
		
		j = First[SAP];
		while(j!=0){
		    int v = Target[j];
		    int future_next  = Next[j];
		    if(v > Nodes_before){
		      j = future_next;
		      continue;
		    }
		    if ( ( r_dlabel2pre[v] >= r_dlabel2pre[SAP] ) && ( r_dlabel2pre[v] < r_dlabel2pre[SAP] + r_dsize[SAP] ) ) {
		      //split
		      delete_arc_from_graph(SAP,v,j);
		      if(r_idom[v] == SAP || Nearest_SAP_Reverse[v] == SAP){
			    insert_arc_to_graph(x,v);
		      }
		      else{
			stack_Cond[v] = true;
		      }
		    }
		    j=future_next;
		}
	    }
	  }
	}
    }
    else{
      int second_vertex;
      if(temp_stack[0] == vertex) second_vertex = temp_stack[1];
      else 			  second_vertex = temp_stack[0];
      
      int Nodes_before = Nodes + j_global;
      testsc(vertex, second_vertex, size_of_SCC);
      
      if(Nodes_before != Nodes + j_global) num_SAP = 1;
    }
    
    return num_SAP;
}

/*breadth-first search from vertex k for delta+1 edges 
   stores parents and labels   */

void BFS (int u, int *Target, int *Next, int *First, bool *marked_S2_Cond){
  
//     printf("========================= (vertices %d / %d)\n",Nodes+j_global, n);
    
    queue_BFS[n]=u;
    marked_BFS[u] = 1;
    int q_out = n;
    int q_in = n-1;
    
    x_sp = u;
    
    while(q_out!=q_in) {
	if (w[queue_BFS[q_out]]> delta && w[parent[queue_BFS[q_out]]] > delta && F_two <= (delta) && queue_BFS[q_out]!= u) {
	      F_two=F_two+1;
	      int v=parent[queue_BFS[q_out]];
	      if (marked_BFS[v] == 0) {
		  marked_BFS[v] = 1;
		  
		  marked_S2[p_S2] = v;
		  marked_S2_Cond[v] = true;
		  p_S2++;
		  
		  queue_BFS[q_in] = v;
		  q_in--;
	      }
	}
	
	int l=First[queue_BFS[q_out]];
	while(l!=0) {
	    int v=Target[l];
	    int future_next = Next[l];
	    
	    if (F_two >(delta)){
		  break;
	    }
	    
	    F_two=F_two + 1;
	    
	    if(marked_BFS[v] == 0 && w[v] > delta && v != u){
		int q = parent[v];
		if (marked_BFS[q] == 0) {
		    marked_BFS[q] = 1;
		    
		    marked_S2[p_S2] = q;
		    marked_S2_Cond[q] = true;
		    p_S2++;
		    
		    queue_BFS[q_in] = q;
		    q_in--;
		    
		}
		if(label2pre[v] > label2pre[x_sp]) x_sp = v;
	    }
	    else if(marked_BFS[v] == 0){
		  marked_BFS[v] = 1;
		  
		  marked_S2[p_S2] = v;
		  marked_S2_Cond[v] = true;
		  p_S2++;
		  
		  queue_BFS[q_in] = v;
		  q_in--;
	      
	    }
	    l = future_next;
	}
	q_out--;
  }
  
  if(!marked_S2_Cond[x_sp]){
	marked_S2[p_S2] = x_sp;
	marked_S2_Cond[x_sp] = true;
	p_S2++;
  }
}




/* depth-first search from vertex k for 2*delta+1 edges 
   stores parents and labels        */

void DFS_1(int u, int *Target, int *Next, int *First, bool *marked_S1_Cond) {
    label2pre[u] = (++next);
    w[u] = 0;
   
    int k=First[u];
    while(k!=0){
	int v=Target[k];
	int future_next = Next[k];
     
	if (F_one >=(check)){
	    break;
	}
	
	F_one=F_one+1;
	w[u]=w[u]+1;
     
	if (!label2pre[v]) { 
	      marked_S1[p_S1]=v;
	      marked_S1_Cond[v]= true;
	      p_S1++;
	      DFS_1(v, Target, Next, First, marked_S1_Cond);
	      w[u] += w[v];
	      parent[v] = u;
	}
        k=future_next;
    }  
}




int outgoing_edge_1;
int outgoing_edge_2;

int OneVertexOut(int u, int *Target, int *Next, int *First, bool *marked_S1_Cond, bool *marked_S2_Cond){
      
    next=0;
    F_one=0;
    marked_S1[1]=u;
    marked_S1_Cond[u] = true;
    p_S1=2;
    parent[u]=0;

    
    DFS_1(u, Target, Next, First, marked_S1_Cond);
 
    if (F_one < check){
      number_of_nodes_of_1egeOut=(p_S1);
      p=p_S1;
      idf=1;
      return idf;
    }
    else{
      
      
// 	int k=0;
	// compute heavy path and store the id of heavy vertices
     
	next_2 = 0;
	F_two = 0;
	p_S2 = 2;
	marked_S2[1] = u;
	marked_S2_Cond[u] = true;
	
	BFS(u, Target, Next, First, marked_S2_Cond);
	
	if (F_two  <= (delta)){
	  number_of_nodes_of_1egeOut=p_S2; 
	  p=p_S2;
	  idf=2;
	   
	  return idf;
	}
	else{   
	  number_of_nodes_of_1egeOut=0;
	  idf = 0;
	  return 0;
	}
    }
} 

int max_stack = 0;

void Process_OneVertexOut(int vertex, int size_of_SCC, int check,int recursionDepth) {

  
  if(recursionDepth > maxRecursionDepth) maxRecursionDepth = recursionDepth;
  numOfRecursionCalls ++;
    
//    printf("Process_OneEdgeOut(%d,%d,%d, %d)\n",vertex,size_of_SCC,check,recursionDepth);
   int component_type;
   int component_type_r;
   
   int v,l,k;
   
  fillReachableStack(vertex, temp_stack);
  
  int Nodes_before = Nodes + j_global;
  
  int num_of_SAPs = find_SAPS(vertex, size_of_SCC);
  //printf(" 0 \n");
  
  if(num_of_SAPs == 0){
        if(temp_stack_top > max2VCSSize) max2VCSSize = temp_stack_top;
        if(temp_stack_top < min2VCSSize) min2VCSSize = temp_stack_top;
        sum2VCSSize += temp_stack_top;
        numOf2VCS++; // number of 2ECC
        return;
  }
  for(int i = Nodes_before+1; i <= Nodes+j_global; i++){
    temp_stack[temp_stack_top++] = i;
  }
  
  
  initReachableToInt(label2pre, 0);
  initReachableToInt(label2pre_2, 0);
  initReachableToInt(marked_S1, 0);
  initReachableToInt(marked_S2, 0);
  initReachableToInt(w, 0);
  initReachableToInt(marked_BFS, 0);
  initReachableToFalse(marked_Cond);
  initReachableToFalse(marked_S1_Cond);
  initReachableToFalse(marked_S2_Cond);

  stack_top = 0;
  for(int i=0 ; i < temp_stack_top; i++){
     if(stack_Cond[temp_stack[i]]) stack[++stack_top] = temp_stack[i];
  }
  
   Nodes_before = Nodes + j_global;
   
   while (stack_top > 0){
	if(stack_top > max_stack) max_stack = stack_top;
	int u = stack[stack_top--];
	if (stack_Cond[u]==true && !in_small_comp[u]){
		  stack_Cond[u]=false;
		  //printf("\n check for component vertex %d \n",u);
		  number_of_nodes_of_1egeOut=0;
		  p_S1 = p_S2 = 0;
		  
		  component_type = OneVertexOut(u, Target, Next, First, marked_S1_Cond, marked_S2_Cond);
		 

		    if (component_type != 0){
		          number_of_components = number_of_components+1;
			  
			  if(number_of_nodes_of_1egeOut > max1VOutSize) max1VOutSize = number_of_nodes_of_1egeOut;
			  if(number_of_nodes_of_1egeOut < min1VOutSize) min1VOutSize = number_of_nodes_of_1egeOut;
			  sum1VOutSize += number_of_nodes_of_1egeOut;
			  numOf1VOut++; // number of 2ECC  
// 			  if(number_of_nodes_of_1egeOut) printf("\n number_of_nodes_of_1egeOut= %d, check=%d \n ",number_of_nodes_of_1egeOut,check);
			  
			  
			  int *S_one;
			  if(idf == 2){
			      S_one= marked_S2;
			  }
			  else{
			      S_one= marked_S1;
			  }
			  
			 //split 
			 
			  int x = 0; 
			  if (idf==2){
				j_global = j_global + 1;
				
				x = Nodes + j_global; 
				real_arc[x] = real_arc[x_sp];
				
				temp_stack[temp_stack_top++] = x;
// 				printf("Added vertex %d/%d\n",temp_stack_top,n);
				
				if(!in_small_comp[x_sp] && stack_Cond[x_sp] == false){
				    //add their endpoints in stack
				    stack[++stack_top]=x_sp;
				    stack_Cond[x_sp]=true;
				}
			  }
			  
			   
			  /// delete inicident edges from the graph 
			  for (int i=1;i<number_of_nodes_of_1egeOut;i++){
				if(idf == 2 && S_one[i] == x_sp) continue;
				l=r_First[S_one[i]];
				while(l!=0){
				      v = r_Target[l]; 
				      int future_next = r_Next[l];
				      
				      if (idf==1){
					if(marked_S1_Cond[v]==false){
					  //printf("\n the edge is from outside the component\n");
					  delete_arc_from_graph(v,S_one[i],l);
					  if(!in_small_comp[v] && stack_Cond[v] == false){
					      //add their endpoints in stack
					      stack[++stack_top]=v;
					      stack_Cond[v]=true;
					  }
					}
				      }
				      else if (idf==2){
					if (v == x_sp){
					    delete_arc_from_graph(v, S_one[i],l);
					    insert_arc_to_graph(x, S_one[i]);
					}   
					if(marked_S2_Cond[v]==false && v != (Nodes + j_global)){
					    delete_arc_from_graph(v,S_one[i],l);
					    if(!in_small_comp[v] && stack_Cond[v] == false){
					      //add their endpoints in stack
					      stack[++stack_top]=v;
					      stack_Cond[v]=true;
					    }
					}
				      }
				      l=future_next;
				}
				if (idf==2){
				    l=First[S_one[i]];
				    while(l!=0){
					  v = Target[l];
					  int future_next = Next[l];
					  if (v == x_sp){
					      delete_arc_from_graph(S_one[i],v,l);
					      insert_arc_to_graph(S_one[i],x);
					  }   
					  l=future_next;
				    }
				}
				stack_Cond[S_one[i]] = false;
				in_small_comp[S_one[i]] = true;
			}
		      
		    }

		    //// initialize the structures that used in OneEdgeOut
		    x_sp = 0;
		    if (p_S1>p_S2){
			k=p_S1;
		    }
		    else{
			k=p_S2;
		    }
		    
		  
		    for (int i=1;i<=k;i++){
		  
			marked_S1_Cond[marked_S1[i]]=false;
			marked_S2_Cond[marked_S2[i]]=false;
			
			marked_BFS [marked_S1[i]] = 0;
			marked_BFS [marked_S2[i]] = 0;
			
			label2pre[marked_S1[i]]=0;
			label2pre_2[marked_S2[i]]=0;
			
			w[marked_S1[i]]=0;
			
			marked_S1[i]=0;
			marked_S2[i]=0;
			
			
			w_matrix[i]=0;
			
		    }
		  
		    if(component_type == 0){
			//////////////////////////////////////////////  
			number_of_nodes_of_1egeOut=0;
			//////////////  REVERSE GRAPH //////////////////////////
			
			p_S1 = p_S2 = 0;
			component_type_r = OneVertexOut(u, r_Target, r_Next, r_First, marked_S1_Cond, marked_S2_Cond);

			if (component_type_r != 0){ 
				  if(number_of_nodes_of_1egeOut > max1VOutSize) max1VOutSize = number_of_nodes_of_1egeOut;
				  if(number_of_nodes_of_1egeOut < min1VOutSize) min1VOutSize = number_of_nodes_of_1egeOut;
				  sum1VOutSize += number_of_nodes_of_1egeOut;
				  numOf1VOut++; // number of 2ECC  
				  number_of_components=number_of_components+1;
				  //if(number_of_nodes_of_1egeOut> 40) printf("\n number_of_nodes_of_1egeOut= %d, check=%d \n ",number_of_nodes_of_1egeOut,check);
				  /// delete inicident edges from the graph 
				  
				  int *S_one_r;
				  
				  if(idf == 2){
				    
				      S_one_r = marked_S2;
				  }
				  else{
				    S_one_r = marked_S1;
				  }
				  
				    int x = 0; 
    				    if (idf==2){
    				      
    				      j_global = j_global + 1;
    				      
    				      x = Nodes + j_global; 
    				      real_arc[x] = real_arc[x_sp];
				      
				      temp_stack[temp_stack_top++] = x;
				      
				      if(!in_small_comp[x_sp] && stack_Cond[x_sp] == false){
					  //add their endpoints in stack
					  stack[++stack_top]=x_sp;
					  stack_Cond[x_sp]=true;
				      }
    				      
    				    }
				  
				  //else{continue;}
				  for (int i=1;i<number_of_nodes_of_1egeOut;i++){
					if(idf == 2 && S_one_r[i] == x_sp) continue;
					l=First[S_one_r[i]];
					while(l!=0){
					      v = Target[l]; 
					      int future_next = Next[l];
					      
					      if (idf==1){
						  if(marked_S1_Cond[v]==false){
						    delete_arc_from_graph(S_one_r[i], v,l);
						    if(!in_small_comp[v] && stack_Cond[v] == false){
							//add their endpoints in stack
							stack[++stack_top]=v;
							stack_Cond[v]=true;
						    }
						    //printf("\n the edge belongs to the component\n");
						    l=future_next;
						    continue;
						  }
					      }
					      else if (idf==2){
						  if (v == x_sp){
						      delete_arc_from_graph(S_one_r[i], v, l);
						      insert_arc_to_graph(S_one_r[i],x);
						  }   
						  if (marked_S2_Cond[v] == false && v != (Nodes + j_global)){
						      delete_arc_from_graph(S_one_r[i],v,l);
						      if(!in_small_comp[v]){
							//add their endpoints in stack
							if(stack_Cond[v] == false){
							  stack[++stack_top]=v;
							  stack_Cond[v]=true;
							}
						      }
						  }   
					      }
					      l=future_next;
					}
					
					if (idf==2){
					    l=r_First[S_one_r[i]];
					    while(l!=0){
						  v = r_Target[l];
						  int future_next = r_Next[l];
						  if (v == x_sp){
						      delete_arc_from_graph(v, S_one_r[i],l);
						      insert_arc_to_graph(x, S_one_r[i]);
						  }   
						  l=future_next;
					    }
					}
					
					stack_Cond[S_one_r[i]] = false;
					in_small_comp[S_one_r[i]] = true;
				    
				  }
			    }  

			  
			  
			  
			  //// initialize the structures that used in OneEdgeOut
			  x_sp = 0;
			  if (p_S1>p_S2){
			      k=p_S1;
			  }
			  else{
			      k=p_S2;
			  }
			  
			
			  for (int i=1;i<=k;i++){
			      marked_S1_Cond[marked_S1[i]]=false;
			      marked_S2_Cond[marked_S2[i]]=false;
			      
			      marked_BFS [marked_S1[i]] = 0;
			      marked_BFS [marked_S2[i]] = 0;
			      
			      label2pre[marked_S1[i]]=0;
			      label2pre_2[marked_S2[i]]=0;
				
			      w[marked_S1[i]]=0;
			      
			      marked_S1[i]=0;
			      marked_S2[i]=0;
			  }
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
	if((SCCs_first[i+1]-SCCs_first[i]) < check){
	    for(int j = SCCs_first[i]; j < SCCs_first[i+1]; j++){ // for every vertex in the SCC
		x = SCCs[j];
		in_small_comp[x] = true;
	    }
	}
    }
    
    for(int i = 1; i < currentSccID; i++){
	for(int j = SCCs_first[i]; j < SCCs_first[i+1]; j++){ // for every vertex in the SCC
	    x = SCCs[j];
	    
	    l=First[x];
	    while(l!=0){
		  v=Target[l]; 
		  int future_next = Next[l];
		  if(SCC_ID[v] != i){
		    delete_arc_from_graph(x,v,l);
		    if(!in_small_comp[v]){
			stack_Cond[v]=true;
		    }
		    if(!in_small_comp[x]){
			stack_Cond[x]=true;
		    }
		  }
		  l=future_next;
	    }
	}
    }

    int *swap_reachable  = new int [temp_stack_top];
    int swap_reachable_top = temp_stack_top;
    Nodes_before = Nodes + j_global;
    
    for(int i = 0; i < temp_stack_top; i++) swap_reachable[i] = temp_stack[i];

    // store the size and a vertex for each SCC (because the SCC matrices are global)
    for(int i = 1; i < currentSccID; i++){
	num_in_SCC[i] = SCCs_first[i+1] - SCCs_first[i];
	vertex_in_SCC[i] = SCCs[SCCs_first[i]];
    }

  
    for(int i = 1; i < currentSccID; i++){
// 	printf("SCC %d / %d with %d vertices, one vertex in the SCC is %d\n",i,currentSccID-1,num_in_SCC[i],vertex_in_SCC[i]);
	if(num_in_SCC[i] > 2){
	    // recursive call
	    (void) find_SAPS(vertex_in_SCC[i], num_in_SCC[i]);
	}
    }

    for(int i = 0; i < swap_reachable_top; i++) temp_stack[i] = swap_reachable[i];
    
    temp_stack_top = swap_reachable_top;

    for(int i = Nodes_before+1; i <= Nodes+j_global; i++){
      temp_stack[temp_stack_top++] = i;
    }
  
    delete [] swap_reachable;
    delete [] num_in_SCC;
    delete [] vertex_in_SCC;
    
    currentSccID = 1;
    connectedComponents(&currentSccID);

    num_in_SCC = new int [currentSccID];
    vertex_in_SCC = new int [currentSccID];

    for(int i = 1; i < currentSccID; i++){
	if((SCCs_first[i+1]-SCCs_first[i]) < check){
	    for(int j = SCCs_first[i]; j < SCCs_first[i+1]; j++){ // for every vertex in the SCC
		x = SCCs[j];
		in_small_comp[x] = true;
	    }
	}
    }
    
    for(int i = 1; i < currentSccID; i++){
	for(int j = SCCs_first[i]; j < SCCs_first[i+1]; j++){ // for every vertex in the SCC
	    x = SCCs[j];
	    l=First[x];
	    while(l!=0){
		  v=Target[l];
		  int future_next = Next[l];
		  if(SCC_ID[v] != i){
		    delete_arc_from_graph(x,v,l);
		    if(!in_small_comp[v]){
			stack_Cond[v]=true;
		    }
		    if(!in_small_comp[x]){
			stack_Cond[x]=true;
		    }
		  }
		  l=future_next;
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
//  	    printf("Call %d-th SCC (out of %d)\n",i,currentSccID-1);
	    Process_OneVertexOut(vertex_in_SCC[i], num_in_SCC[i], check,recursionDepth+1);
	}
    }
    
    delete [] num_in_SCC;
    delete [] vertex_in_SCC;
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
     delete [] label2pre_2; 
     delete [] marked_Cond ;
   
     delete [] marked_S1_Cond;
     delete [] marked_S2_Cond;
    
     delete [] stack_Cond;
     delete [] in_small_comp;
    
     delete [] w ;
    
     delete [] LowLink ;
    
     delete [] marked_S1;
     delete [] marked_S2;
     delete [] stack ;
     delete [] temp_stack ;
     delete [] temp_inserted ;
    
     delete [] SCCs;
     delete [] SCCs_first;
     delete [] SCC_ID ;
     
     delete [] queue_BFS;
     delete [] w_matrix;
     delete [] real_arc;
}





void init() {
  
    n=2*Nodes;
    m=M;
    
    
    First=new int [n+2];
    Target=new int [m+2];
    Next=new int [m+2];
    Before=new int [m+2];
    next_free = 1;
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
    
    LowLink = new int [n+2];
    queue_BFS = new int [n+2];
    
    marked_S1=new int [n+2];
    marked_S2=new int [n+2];
    stack = new int [n+2];
    temp_stack = new int [n+2];
    temp_inserted = new bool[n+2];
    marked_BFS = new int [n+2];
    queue_BFS = new int  [n+2];
    
    SCCs = new int [n+2];
    SCCs_first= new int [n+2];
    SCC_ID = new int [n+2];
    w_matrix = new int [n+2]; 
    real_arc= new int [n+2];
    
    SAP_in_Forward = new bool [n+2];
    SAP_in_Reverse = new bool [n+2];
    
    max2VCSSize = sum2VCSSize = numOf2VCS = maxRecursionDepth = numOfRecursionCalls = 0;
    max1VOutSize=sum1VOutSize = numOf1VOut = 0;
    min2VCSSize=min1VOutSize=n;

    for (int i=0;i<n+2;i++){
      
	stack[i]=0;
	queue_BFS[i]=0;
	
	First[i]=0; 
	r_First[i]=0;
	
	parent[i] = 0;
	
	in_small_comp[i] = false;
	stack_Cond[i]=true;
	temp_inserted[i] = false;
	marked_BFS [i] = 0;
	marked_S1 [i] = 0;
	marked_S2 [i] = 0;
	marked_S1_Cond[i] = false;
	marked_S2_Cond[i] = false;
	w_matrix[i] = 0;
	real_arc[i] = i;	
	
    }
   
    
    for (int i=0;i<m+2;i++){
      Before[i]=0;
      r_Next[i]=0;
      r_Before[i]=0;
    }
    for (int i=1;i<m+1;i++){
      Next[i]=i+1;
    }
    Next[0] = Next[M+1] = 0;

  
}



void builtgraph() {
    int input_source,input_target;
    
    for (int current_input_pos=0 ; current_input_pos < M; current_input_pos++){
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
            case 'p': if (sscanf(line, "p %d %d", &Nodes, &M) != 2) {
                    fprintf(stderr, "Error reading graph size (%s).\n", file);
                    exit(-1);
                }
                input_file = new int [2*M];
                break;
            case 'a': if (sscanf(line, "a %d %d", &x, &y) != 2) {
                    fprintf(stderr, "Error reading graph arc (%s).\n", file);
                    exit(-1);
                }
                assert(x <= Nodes);
                assert(y <= Nodes);

                if(x == y) break;

                input_file[p++] = x;
                input_file[p++] = y;
                if (p>2*M) {
                    fprintf(stderr, "Error! Graph has >%d arcs.\n", M);
                    exit(-1);
                }
                break;
            default: fprintf(stderr, "Error reading graph (%s).\n", file);
                exit(-1);
        }
    }
    printf("number of nodes =%d , number of edges=%d\n",Nodes,M);
    fprintf(stderr, "END reading graph (%s).\n", file);
    fclose(input);
}



int main(int argc, char *argv[]) {
    if (argc != 4 ) {
        printf("\n usage: %s <input file> <seed> <starting vertex> <parameter for component> \n\n", argv[0]);
        exit(-1);
    }

    char *file = argv[1];
    int seed = atoi(argv[2]);
    int rand_node = atoi(argv[3]);
    

    readGraph(file);
  
    //To Create the output File
    inputFileName = Utilities::TO_STRING(file);
    preText = Utilities::GetFileName(inputFileName);
    testDate = Utilities::GetCurrentDateTime();
    
    RFWTimer timer(true);
    
    init();
    
    builtgraph();
     
  //  output = fopen("exit.txt", "w");
    
    srand(seed);
    delta=sqrt(M); 
    check=2*delta + 1; // check=2*delta +1 
    //slt(1,N);
    k_split = 0;
    Process_OneVertexOut(rand_node,Nodes,check,1);
    delete_space();
    
    takenTime = timer.getTime();
   
    ShowOutput();
    if (is_export_output)
        ExportOutput();
    
//     //if(numOf2ECS > 0){
//         fprintf (stdout, "Number of 2VCCs       : %d\n", numOf2VCS);
//         fprintf (stdout, "Average 2VCC size     : %f\n", (float)sum2VCSSize / (float)numOf2VCS);
//         fprintf (stdout, "Max 2VCC              : %d\n", max2VCSSize);
//         fprintf (stdout, "Min 2VCC              : %d\n", min2VCSSize);
//    // }
//     
//     fprintf (stdout,     "====================================\n");
//     fprintf (stdout,     "Total calls to main		: %d\n", numOfRecursionCalls);
//     fprintf (stdout,     "Max recursion depth		: %d\n", maxRecursionDepth);
//     
//     fprintf (stdout,     "====================================\n");
//     fprintf (stdout, "Delta = %d, 2*Delta+1 = %d\n", delta, check);
//     fprintf (stdout, "Number of 1EOut components       : %d\n", numOf1VOut);
//     fprintf (stdout, "Average 1EOut component size     : %f\n", (float)sum1VOutSize / (float)numOf1VOut);
//     fprintf (stdout, "Max 1VOut component              : %d\n", max1VOutSize);
//     fprintf (stdout, "Min 1VOut component              : %d\n", min1VOutSize);
//     
//     fprintf (stdout, "totaltime %f\n", t);
    
}



// int maxRecursionDepth, max2VCSSize, min2VCSSize, sum2VCSSize, numOf2VCS, numOfRecursionCalls;
// int max1VOutSize, min1VOutSize, sum1VOutSize, numOf1VOut;
void ShowOutput() {

    // Source code added by Nilu

    cout << " x==========================================x" << endl;
    cout << " |     Summary of the 2VSCC Computation     |" << endl;
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
    cout << " Total 2VSCC components   : " << numOf2VCS << endl;
    cout << " Total Time Taken         : " << takenTime << " seconds" << endl;
    cout << " Recursion Depth Level    : " << maxRecursionDepth << endl;
    cout << " Minimum 2VSCC Vertex Size: " << min2VCSSize << endl;
    cout << " Maximum 2VSCC Vertex Size: " << max2VCSSize << endl;
    cout << " Average 2VSCC Vertex Size: " << avg2eccSize << endl;
    cout << " Total Fractions          : " << numOfRecursionCalls << endl;

}

void ExportOutput() {

    string fileName = preText + "_2VSCC_TIME_" + algorithm_code + ".txt";
    ofstream out_file;
    out_file.open(fileName.c_str());

    //Header on file
    out_file << " x==========================================x" << endl;
    out_file << " |     Summary of the 2VSCC Computation     |" << endl;
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
    out_file << " Total 2VSCC components   : " << numOf2VCS << endl;
    out_file << " Total Time Taken         : " << takenTime << " seconds"
            << endl;
    out_file << " Total Recursion Depth    : " << maxRecursionDepth << endl;

    out_file << " Minimum 2VSCC Vertex Size: " << min2VCSSize << endl;
    out_file << " Maximum 2VSCC Vertex Size: " << max2VCSSize << endl;
    out_file << " Average 2VSCC Vertex Size: " << avg2eccSize << endl;
    out_file << " Total Fractions          : " << numOfRecursionCalls << endl;

    out_file.close();

}
