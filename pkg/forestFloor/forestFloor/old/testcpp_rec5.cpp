#include <Rcpp.h>
using namespace Rcpp;

//defining recursive function to iterate nodes
void follow_path(
//local recursive function environment                 
                 bool calculate_node_pred,
                 int i_tree,
                 int this_node,
                 double parent_pred,
                 int parent_bestvar,
                 int passed_OOB_count,
                 int passed_IB_count,
                 IntegerVector passed_innodes,
                 IntegerVector train_innodes,
//global R-objects
                 NumericMatrix X,                  // 1  X
                 NumericVector Y,
                 IntegerMatrix leftDaughter,       // 6  LD
                 IntegerMatrix rightDaughter,      // 7  RD
                 IntegerMatrix nodestatus,         // 8  nodestatus
                 NumericMatrix xbestsplit,         // 10 xsplits
                 NumericMatrix nodepred,           // 11 averagetrainnodes
                 IntegerMatrix bestvar,            //12 bestvar
                 NumericMatrix localIncrements) 
{

 //, IntegerVector parent_innode, int parent_nOOBs) {
//computing and adding the latest increment to matrix
double current_pred = 0;
//printf("passed_B_count %d \n",passed_IB_count);
if(calculate_node_pred) {
  for(int thisIB_count = 0; thisIB_count<passed_IB_count;thisIB_count++) {
    current_pred += Y(train_innodes(thisIB_count));
    //printf("thisIB_ind %d \n",train_innodes(thisIB_count));
  }
  current_pred /= (passed_IB_count); // no +1 to denominator
  //printf("current_pred %f", current_pred);
  
} else { 
  current_pred = nodepred(this_node,i_tree); //reuse node_pred from RF-object, only regression
}
  double this_increment = current_pred - parent_pred;
  int outcome = 0;
  int this_obs=0;
  



for(int i_obs=0;i_obs<passed_OOB_count;i_obs++) {
    this_obs = passed_innodes[i_obs];
    localIncrements(this_obs,parent_bestvar) += this_increment;  
}

//Rprintf("in node %d increment %d \n",this_node,this_increment);

//localIncrements(1,1) = 49;
//if(this_node==-61) return 57;
  
  
  //if not a terminal node, split 
  if(passed_OOB_count>0) {
    if(nodestatus(this_node,i_tree)==-3) { // #if this is not a terminal node
      int current_bestvar = bestvar(this_node,i_tree);
      int OOB_count_left  = 0;
      int OOB_count_right = 0;
      int IB_count_left   = 0;
      int IB_count_right  = 0;
      IntegerVector OOBs_leftnode (passed_OOB_count);
      IntegerVector OOBs_rightnode(passed_OOB_count);
      IntegerVector IBs_leftnode (passed_IB_count);
      IntegerVector IBs_rightnode(passed_IB_count);
      
      double this_split = xbestsplit(this_node,i_tree);
      //splitting OOB
      for(int i_obs = 0;i_obs<passed_OOB_count;i_obs++) {
        if(X(passed_innodes[i_obs],current_bestvar) < this_split) {
          OOBs_leftnode[OOB_count_left] = passed_innodes[i_obs];
          OOB_count_left++;
        } else {
          OOBs_rightnode[OOB_count_right] = passed_innodes[i_obs];
          OOB_count_right++;
        }
      }
      
      //splitting IB
      //printf("passed ib count is %d \n", passed_IB_count);
      for(int i_obs = 0;i_obs<passed_IB_count;i_obs++) {
        if(X(train_innodes[i_obs],current_bestvar) < this_split) {
          IBs_leftnode[IB_count_left] = train_innodes[i_obs];
          IB_count_left++;
        } else {
          IBs_rightnode[IB_count_right] = train_innodes[i_obs];
          IB_count_right++;
        }
      }
      //printf("got to here, left %d right %d \n",IB_count_left,IB_count_right);
      //for(int i_print=0;i_print<IB_count_left;i_print++) printf("%d %d \n ",i_print,IBs_leftnode(i_print));
      
      
      //initiate left step
      if(OOB_count_left>0) {
      follow_path(
//local recursive function environment                 
        calculate_node_pred,
        i_tree,
        leftDaughter(this_node,i_tree),
        current_pred,
        current_bestvar,
        OOB_count_left,
        IB_count_left,
        OOBs_leftnode,
        IBs_leftnode,
//pointers to global R-objects
        X,                  // 1  X
        Y,
        leftDaughter,       // 6  LD
        rightDaughter,      // 7  RD
        nodestatus,         // 8  nodestatus
        xbestsplit,         // 10 xsplits
        nodepred,           // 11 averagetrainnodes
        bestvar,            //12 bestvar
        localIncrements);
        OOB_count_left=0;
        //Rprintf("going back from left, now in %d \n",this_node);
        //if(outcome != 0) return outcome;
      }
      
      //... came back from a left step, now going right
      if(OOB_count_right>0) {
        follow_path(
        calculate_node_pred,
        i_tree,
        rightDaughter(this_node,i_tree),
        current_pred,
        current_bestvar,
        OOB_count_right,
        IB_count_right,
        OOBs_rightnode,
        IBs_rightnode,
        X,                  // 1  X
        Y,
        leftDaughter,       // 6  LD
        rightDaughter,      // 7  RD
        nodestatus,         // 8  nodestatus
        xbestsplit,         // 10 xsplits
        nodepred,           // 11 averagetrainnodes
        bestvar,            //12 bestvar
        localIncrements);
        OOB_count_right=0;
        //Rprintf("going back from righ, now in %d \n",this_node);
        //if(outcome != 0) return outcome;
        }
        
        
      //... came back from a right step, all nodes below have been completed
        nodestatus(this_node,i_tree)==-1;
    }
  }
  //got to here as node was terminal or all nodes below have been checked
  //Now leaving node going up. Next OOB indices is saved locally in upper node as OOBs_leftnode or OOBs_rightnode
  // if upper node is rootnode, the path will be terminated
  //no returns all changes were applied to global NumericVector localIncrements
//return 0;
}


/// defining Rcpp function to communicate with R
//[[Rcpp::export]]
void recTree(int  vars,               //local 3  nvar
            int  obs,                 //local 4  nobs
            int  ntree,               //local  5  ntrees
            bool calculate_node_pred, //should node prediction
            NumericMatrix X,                  // 1  X
            NumericVector Y,
            IntegerMatrix leftDaughter,       // 6  LD
            IntegerMatrix rightDaughter,      // 7  RD
            IntegerMatrix nodestatus,         // 8  nodestatus
            NumericMatrix xbestsplit,         // 10 xsplits
            NumericMatrix nodepred,           // 11 averagetrainnodes
            IntegerMatrix bestvar,            // 12 bestvar
            IntegerMatrix inbag,
            NumericMatrix localIncrements)    // 15 inbag obsXtrees
{
  
  
  //declare internal function variables
  
  IntegerVector innodes_root(obs);        //list of OOBs 
  IntegerVector train_innodes_root(obs);  //list of IBs
  int this_bagcount =0;
  int times_inbag = 0;
  int OOB_count = 0;
  int IB_count = 0;
  int outcome = 0;
  double root_pred = 0;
  
  //iterate each tree and compute and sum to localIncrements
  for(int i_tree=0;i_tree<ntree;i_tree++){
    //make range of Out Of Bag observations in root of tree
    OOB_count = 0;//reset OOB_count for this new tree
    IB_count =0 ; 
    root_pred =  0;
    for(int i_obs=0;i_obs<obs;i_obs++) { 
      this_bagcount=inbag(i_obs,i_tree);
      if(this_bagcount==0) {      // if observation was not inbag..
        innodes_root[OOB_count] = i_obs;  // add observation_indice to list
        OOB_count++;
      } else {
        if(calculate_node_pred) {
          //printf("obs no. %d is was used %d \n",i_obs,this_bagcount);
          for(times_inbag = this_bagcount;times_inbag>0;times_inbag--) {
            //printf("times_inbag %d \n",times_inbag);
            train_innodes_root[IB_count] = i_obs;
            root_pred += Y(i_obs);
            IB_count++;
            //printf("ibcount %d",IB_count);
            if(IB_count>(obs)) {
              //printf("error maxout IB_count");
              return;
            }
          }
          
  
        }  
      }
      
    }
    root_pred /= IB_count;
    if(!calculate_node_pred) root_pred= nodepred(0,i_tree);
    //printf("root pred %f \n",root_pred);

    //initiating varibles for recursive search of tree
  
     follow_path(
//local recursive function environment                 
        calculate_node_pred,
        i_tree,
        0,   //start in root node, 0
        root_pred, //parent_pred set to root_pred,
        0,    //dummy number, any var within X vars will do 
        OOB_count,  //how many obs OOB to start with in root node
        IB_count,
        innodes_root,  //X.rows indices of OOB observations in rootnode
        train_innodes_root,
//pointers to global R-objects
        X,                  // 1  X dataset
        Y,
        leftDaughter,       // 6  LD
        rightDaughter,      // 7  RD
        nodestatus,         // 8  nodestatus
        xbestsplit,         // 10 xsplits
        nodepred,           // 11 averagetrainnodes
        bestvar,            //12 bestvar
        localIncrements);
       // if(outcome != 0) return outcome;
  } // go to next tree
  
  //divide sum of increments with trees iterated
  for(int i_obs=0;i_obs<obs;i_obs++){
    for(int i_vars=0;i_vars<vars;i_vars++){
      localIncrements(i_obs,i_vars) /= ntree;
    }
  }

}