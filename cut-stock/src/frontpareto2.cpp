/*-------+-------------------------------------------------------------------------+------------+
         | See file LICENSE at the root of the git project for licence information |
         +------------------------------------------------------------------------*/

#include<cstdio>
#include<cstring>
#include<cassert>
#include<cstdlib>
#include<iostream>
#include<ctime>
using namespace std;

#include "frontpareto2.h"
#define ROOT_NODE 1
#define LEFT_NODE 2
#define RIGHT_NODE 3


bst* frontpareto2::new_leaf(cst_t c, val_t v, bst* up, int type){
    bst* tmp  = new bst;
    tmp -> cst = c;
    tmp -> val = v;
    tmp -> left   = NULL;
    tmp -> right  = NULL;
    tmp -> up     = up;
    tmp -> type   = type;
    #ifdef ATTACH_INFO_TO_PAIRS
    lastAdded = tmp;
    #endif
    return tmp;
}
bst* frontpareto2::dive_left(bst* r){
    while(r->left!=NULL)
        r = r->left;
    return r;
}

bst* frontpareto2::dive_right(bst* r){
    while(r->right!=NULL)
        r = r->right;
    return r;
}

bst* frontpareto2::next(bst* r){
    if(r->right != NULL)              //if right child exists
        return dive_left(r->right);   //reach the left-most elem on right branch
    if(r->type == LEFT_NODE)          //if no right child but we are on a left node 
        return r->up;                 //we can go up and that is sufficient
    //We get here: we are on right branch
    //Move up until we get into a node that is a left branch or root
    while(r->type == RIGHT_NODE){
        r = r-> up;
    }
    return r-> up;
}

bst* frontpareto2::prev(bst* r){
    if(r->left != NULL)               //if left child exists
        return dive_right(r->left);   //reach the right-most elem on left branch
    if(r->type == RIGHT_NODE)         //if no left child but we are on a right node 
        return r->up;                 //we can go up and that is sufficient
    //We get here: we are on left branch
    //Move up until we get into a node that is a right branch or root
    while(r->type == LEFT_NODE){
        r = r-> up;
    }
    return r-> up;
}


void frontpareto2::remove_node(bst* r){
    if(r->right == NULL){
        if(r->type==ROOT_NODE){
            r = r-> left;
            delete(root);
            root = r;
            root->type = ROOT_NODE;
            root->up   = NULL;
            return;
        }
        bst* parent  = r->up;
        assert(parent!=NULL);
        int type     = r->type;
        bst* child   = r->left;
        if(child != NULL){//r is not a leaf
            child-> up   = parent;
            child-> type = type;
        }
        delete (r);
        if(type==LEFT_NODE)
            parent->left = child;
        else
            parent->right = child;
    }else{
        int type        = r->type;
        bst* parent     = r->up;
        bst* left       = r->left;
        bst* right      = r->right;
        bst* rightleft  = dive_left(r->right);
        rightleft->left = left;
        if(r->left!=NULL)
            r->left->up = rightleft;
        if(r->right!=NULL)
            r->right->up = rightleft;
        delete (r);
        r = right;
        r->type = type;
        if(type==ROOT_NODE){
            root = r;
            root-> up = NULL;
        }
        else{
            r->up = parent;
            if(type==LEFT_NODE)
                parent->left = r;
            else
                parent->right = r;
        }
    }
}
//remove all states that become dominated after inserting r
void frontpareto2::filter(bst* r){
    val_t val_new = r->val;
    bst* iter = next(r);
    while( (iter!=NULL) && (iter->val <= val_new)) {
        bst* iter2 = next(iter);
        remove_node(iter);
        iter = iter2;
    }
}
int frontpareto2::addIfHigherVal(bst* r, cst_t c, val_t v){
    if(r->cst <= c){
        if(r->val >= v)
            return 0;
        if(r->cst==c)
        if(r->val < v){
            r->val = v;
            filter(r);
            return 1;
        }
        if(r->right != NULL) 
            return addIfHigherVal(r->right,c,v);
        else{
            r->right = new_leaf(c,v,r,RIGHT_NODE);
            filter(r->right);
            return 1;
        }
    }else{
        if(r->left != NULL) 
            return addIfHigherVal(r->left,c,v);
        else{
            r->left = new_leaf(c,v,r,LEFT_NODE);
            filter(r->left);
            return 1;
        }
    }
}
int frontpareto2::addIfHigherVal(cst_t c, val_t v){
    if(root==NULL){
        root = new_leaf(c,v,NULL,ROOT_NODE);//NULL=no parent
        return 1;
    }
    return addIfHigherVal(root, c, v);
}

void frontpareto2::freeMem(bst*r){
    if(r->left!=NULL)
        freeMem(r->left);
    if(r->right!=NULL)
        freeMem(r->right);
    delete r;
}
void frontpareto2::freeMem(){
    freeMem(root);
}
void frontpareto2::liste_print(bst* r, char* prefix, int left_child){
    if(r==NULL) 
        return;
    char* new_prefix = new char[strlen(prefix)+5+1];
    sprintf(new_prefix, "%s    |", prefix);
    if(strlen(prefix)>0)
        if(left_child)
            new_prefix[strlen(prefix)-1]=' ';
    liste_print(r->left, new_prefix,1);
    if(strlen(prefix)>0)
        if(left_child)
            new_prefix[strlen(prefix)-1]='|';

    if(strlen(prefix)>0){
        char tmp = prefix[strlen(prefix)-1];
        prefix[strlen(prefix)-1]='+';
        cout<<prefix<<"---"<< r->cst <<"/"<< r->val <<endl;
        prefix[strlen(prefix)-1] = tmp;
    }else
        cout<< r->cst << "/" << r->val <<endl;

    if(strlen(prefix)>0)
        if(!left_child)
            new_prefix[strlen(prefix)-1]=' ';
    liste_print(r->right, new_prefix,0);
    if(strlen(prefix)>0)
        if(!left_child)
            new_prefix[strlen(prefix)-1]='|';
    delete(new_prefix);
    
}

void frontpareto2::print_tree(){
    liste_print(root,(char*)"",0);
}

int frontpareto2::first(cst_t& c, val_t& v){
    if(root==NULL)
        return 0; //Failure
    iterator = dive_left(root);
    c = iterator->cst;
    v = iterator->val;
    return 1;
}
int frontpareto2::next(cst_t& c, val_t& v){
    iterator = next(iterator);
    if(iterator==NULL)
        return 0;
    c = iterator->cst;
    v = iterator->val;
    return 1;
}
int frontpareto2::first(val_t& v){
    if(root==NULL)
        return 0; //Failure
    iterator = dive_left(root);
    v = iterator->val;
    return 1;
}
int frontpareto2::next(val_t& v){
    iterator = next(iterator);
    if(iterator==NULL)
        return 0;
    v = iterator->val;
    return 1;
}
int frontpareto2::last(cst_t& c, val_t& v){
    if(root==NULL)
        return 0; //Failure
    iterator = dive_right(root);
    c = iterator->cst;
    v = iterator->val;
    return 1;
}
int frontpareto2::prev(cst_t& c, val_t& v){
    iterator = prev(iterator);
    if(iterator==NULL)
        return 0;
    c = iterator->cst;
    v = iterator->val;
    return 1;
    
}

#ifdef ATTACH_INFO_TO_PAIRS
void frontpareto2::putInfoOnLastAdded(void* newinfo)
{
    lastAdded->info = newinfo;
}
void* frontpareto2::getInfoCurrElem(){
    return iterator->info;
}

void frontpareto2::freeMem(void(*newCallbackFreeAttached)(void*), bst*r){
    if(r->left!=NULL)
        freeMem(newCallbackFreeAttached,r->left);
    if(r->right!=NULL)
        freeMem(newCallbackFreeAttached,r->right);
    newCallbackFreeAttached(r->info);
    delete r;
}
void frontpareto2::freeMem(void(*newCallbackFreeAttached)(void*)){
    freeMem(newCallbackFreeAttached, root);
}

#endif

//How to use it:
//int main(){
//    frontpareto2 liste;
//    srand(7);
//    cout<<"Here is the tree: "<<endl;
//    for(long i=0;i<91;i++){
//        int p = (int)(rand());
//        double v = (double)((rand()%100));
//        p = p%1000;
//        liste.addIfHigherVal(p,v);
//    }
//    int cost;
//    double v;
//    liste.first(cost,v);
//    liste.print_tree();
//    cout<<"\nThe first is "<<cost<<"/"<<v<<endl;
//    liste.freeMem();
//}
