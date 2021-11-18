/*--------------------------------------------------------------+
| Author: @ Daniel Porumbel 2021                                |
|License: Any person obtaining a copy of this code is free to   |
|         use it in any manner, subject to two conditions:      |
|           1) no profit may ever be made from using this code  |
|           2) these 5 lines of text shall be included          |
+--------------------------------------------------------------*/

#include"frontpareto1.h"
#include<iostream>
#include<set>
#include<climits>
#include<cassert>
using namespace std;

#define DEBCOUT(x) ;
    
//#define DEBCOUT(x) cout<<x;  //use this to see some internal messages
#define COUT(x) cout<<x;

#define MAX_CST_VAL INT_MAX
//The member len means that the interval [cst, cst+len) contains no other element
//The next element could be at point cst+len and has a higher val than this.val

//The use of this option may slow down things a bit, but it can be useful
//if you need attach some additional information to each pareto front element,
//e.g., a precedence relation to reconstruct an optimal solution in the end
#ifndef ATTACH_INFO_TO_PAIRS
struct cell{
    cst_t cst;
    val_t val;
    cst_t len;
    cell(int cst_, val_t val_){
        cst = cst_;
        val = val_;
        len = 0;
    }
};
#else
struct cell{
    cst_t cst;
    val_t val;
    cst_t len;
    void* info;
    cell(cst_t cst_, val_t val_){
        info = NULL;
        cst = cst_;
        val = val_;
        len = 0;
    }
};
struct infoList{
    void* inf;
    infoList* nxt;
};
#endif
struct classcomp {
  //this function emulates an < operator.
  //doc: Two set elements are considered equivalent if < returns false reflexively
  //notice < returns false reflexively if old.cst<=new.cst && old.val>=new.val  (*)
  //                           && there is no elem between old.val and new.val
  //         => a new element is simply rejected when (*) holds.
  bool operator() (cell* const & lhs, cell* const & rhs) const {
        if(lhs->len==0)                             //left is newly added, new elem
            return lhs->cst<rhs->cst;               //is not equiv to elems at right
        if(rhs->len==0) {                           //right is newly added
            if(lhs->cst<=rhs->cst) {
                if(lhs->val>=rhs->val)
                   if(lhs->cst+lhs->len > rhs->cst) //return false for above (*)
                       return false;
                return true;
            }
            return false;                           //right cst < left cst (cst)
        }
        return lhs->cst<rhs->cst;                   //left&right not newly added
  }
  //bool operator() (cell* const & lhs, cell* const & rhs) const{
  //      if(lhs->len==0){
  //          if(rhs->cst<=lhs->cst){
  //          //  if(rhs->val>=lhs->val)
  //          //      if(rhs->cst+rhs->len-1>=lhs->cst)
  //          //          return false;
  //              return false;
  //          }
  //          //if(rhs->cst>lhs->cst)
  //          return true;
  //      }
  //      if(rhs->len==0){
  //          if(lhs->cst<=rhs->cst){
  //              if(lhs->val>=rhs->val)
  //                  if(lhs->cst+lhs->len-1>=rhs->cst)
  //                      return false;
  //              return true;
  //          }
  //          //lhs->cst>rhs->cst
  //          return false;
  //      }
  //      return lhs->cst<rhs->cst;
  //  }
};

class dpointer{
    public:
        //Doc: Sets are typically implemented as binary search trees.
        //Doc: Internally, sets keep all their elements sorted following the criterion specified by its comparison object. The elements are always inserted in its respective position following this ordering.
        set<cell*,classcomp> front;         
        set<cell*,classcomp>::iterator mainIter;
        void(*callbackOnDelete)(cst_t, val_t);
        #ifdef ATTACH_INFO_TO_PAIRS
        cell* lastAdded;
        infoList* deletedInfo;
        infoList* deletedInfoIter;
        void(*callbackFreeAttached)(void*);
        dpointer(){
            deletedInfo          = NULL;
            callbackOnDelete     = NULL;
            callbackFreeAttached = NULL;
        }
        #else
        dpointer(){
            callbackOnDelete = NULL;
        }
        #endif
};


frontpareto1::frontpareto1():dptr(new dpointer()){};
frontpareto1::~frontpareto1(){delete dptr;};
int frontpareto1::addIfHigherVal(cst_t cst, val_t v)
{
    DEBCOUT("********"<<endl);
    set<cell*,classcomp>::iterator it;
    std::pair<std::set<cell*>::iterator,bool> ret;
    cell* theNewCell = new cell(cst,v);
    ret = dptr->front.insert(theNewCell);
    if(ret.second==false){
        //DEBCOUT("doubled:"<<ret.first->cst<<"val="<<ret.first->val<<endl;
        DEBCOUT("could not add "<<cst<<" "<<v<<" because of "<<
            (*ret.first)->cst<<" "<<(*ret.first)->val<<endl);
        delete theNewCell;
        return 0;
    }
    //did find a place for the new element
    DEBCOUT("added"<<(*ret.first)->cst<<","<<(*ret.first)->val);
    it = ret.first;
    #ifdef ATTACH_INFO_TO_PAIRS
    dptr->lastAdded = *it;
    #endif
    //dealing with the element before
    --it;
    if(it!=dptr->front.end()&&(it!=ret.first)){
        (*it)->len= -(*it)->cst + (*ret.first)->cst ;
        DEBCOUT(" after"<<(*it)->cst<<","<<(*it)->val<<endl);
        if((*it)->len==0){
            DEBCOUT("delete existing thing at same position:"<<(*it)->cst<<" "<<(*it)->val<<endl);
            if(dptr->callbackOnDelete!=NULL)
                dptr->callbackOnDelete((*it)->cst,(*it)->val);
            #ifdef ATTACH_INFO_TO_PAIRS
            keepDeleted((*it)->info);
            #endif
            delete *it;
            dptr->front.erase(it);
        }
    }
    else{
        DEBCOUT(" at dptr->front\n");
    }
    it = ret.first;
    it++;
    //dealing with the element after
    set<cell*,classcomp>::iterator it2;
    do{
        if(it==dptr->front.end()){
            (*ret.first)->len=MAX_CST_VAL/2;//redefine MAX_CST_VAL if, eg, cst_t is char
            break;
        }
        if((*it)->val>v){                //next is "ok"
            (*ret.first)->len= (*it)->cst - (*ret.first)->cst ;
            break;
        }
        it2 = it;
        it2++;
        DEBCOUT("and deleted"<<(*it)->cst<<" "<<(*it)->val<<endl);
        if(dptr->callbackOnDelete!=NULL)
            dptr->callbackOnDelete((*it)->cst,(*it)->val);
        #ifdef ATTACH_INFO_TO_PAIRS
        keepDeleted((*it)->info);
        #endif
        delete *it;
        dptr->front.erase(it);
        it  = it2;
    }while(true);
    return 1;
}

void frontpareto1::freeMem()
{
    set<cell*,classcomp>::iterator it;
    for (it=dptr->front.begin(); it!=dptr->front.end(); ++it) {
        #ifdef ATTACH_INFO_TO_PAIRS
        if(((*it)->info)!=NULL){
          DEBCOUT("EErasing info"<<*(int*)((*it)->info));         //needed to cast it to something to print it
          if(dptr->callbackFreeAttached==NULL) delete (*it)->info;//it should never call this, use setCallbackFreeAttachInfo
                      else dptr->callbackFreeAttached((*it)->info);
        }
        #endif
        delete *it;
    }
    dptr->front.clear();
    #ifdef ATTACH_INFO_TO_PAIRS
    while(dptr->deletedInfo!=NULL){
          infoList* toErase = dptr->deletedInfo;
          dptr->deletedInfo = dptr->deletedInfo->nxt;
          DEBCOUT("Erasing info"<<*((int*)toErase->inf));          //needed to cast it to something to print it
          if(dptr->callbackFreeAttached==NULL) delete toErase->inf;//it should never call this, use setCallbackFreeAttachInfo
                      else dptr->callbackFreeAttached(toErase->inf);
          delete toErase;
    }
    #endif
}
void frontpareto1::setDeleteCallback(void(*newCallbackOnDelete)(cst_t, val_t))
{
    dptr->callbackOnDelete = newCallbackOnDelete;
}
#ifdef ATTACH_INFO_TO_PAIRS
void frontpareto1::setCallbackFreeAttachInfo(void(*newCallbackFreeAttached)(void*))
{
    dptr->callbackFreeAttached = newCallbackFreeAttached;
}
#endif
void frontpareto1::printPairs()
{
    set<cell*,classcomp>::iterator it;
    for (it=dptr->front.begin(); it!=dptr->front.end(); ++it) 
        cout<<(*it)->cst << " "<<(*it)->val<< " "<<(*it)->len<<endl;
}
int frontpareto1::first(val_t &vv)
{
    dptr->mainIter = dptr->front.begin();
    if(dptr->mainIter==dptr->front.end())
        return 0;
    vv=(*(dptr->mainIter))->val;
    return 1;
}
int frontpareto1::next(val_t &vv)
{
    ++dptr->mainIter;
    if(dptr->mainIter==dptr->front.end())
        return 0;
    vv=(*(dptr->mainIter))->val;
    return 1;
}
int frontpareto1::first(cst_t &ii, val_t &vv)
{
    dptr->mainIter = dptr->front.begin();
    if(dptr->mainIter==dptr->front.end())
        return 0;
    ii=(*(dptr->mainIter))->cst;
    vv=(*(dptr->mainIter))->val;
    return 1;
}
int frontpareto1::next(cst_t &ii, val_t &vv)
{
    ++dptr->mainIter;
    if(dptr->mainIter==dptr->front.end())
        return 0;
    ii=(*(dptr->mainIter))->cst;
    vv=(*(dptr->mainIter))->val;
    return 1;
}
#ifdef ATTACH_INFO_TO_PAIRS
void frontpareto1::freeMem(void(*newCallbackFreeAttached)(void*))
{
    set<cell*,classcomp>::iterator it;
    for (it=dptr->front.begin(); it!=dptr->front.end(); ++it) {
        if(((*it)->info)!=NULL){
          DEBCOUT("EErrasing info"<<*(int*)((*it)->info));  //needed to cast it to something to print it
          newCallbackFreeAttached((*it)->info);
        }
        delete *it;
    }
    dptr->front.clear();
    while(dptr->deletedInfo!=NULL){
          infoList* toErase = dptr->deletedInfo;
          dptr->deletedInfo = dptr->deletedInfo->nxt;
          DEBCOUT("Errasing info"<<*((int*)toErase->inf)); //needed to cast it to something to print it
          newCallbackFreeAttached(toErase->inf);
          delete toErase;
    }
}
void frontpareto1::keepDeleted(void*currInfo)   
{
   if(currInfo==NULL)
      return;
   if(dptr->deletedInfo==NULL){
       dptr->deletedInfo      = new infoList;
       dptr->deletedInfo->inf = currInfo;
       dptr->deletedInfo->nxt = NULL;
       dptr->deletedInfoIter  = dptr->deletedInfo;
       return;
   }
   assert(dptr->deletedInfoIter->nxt==NULL);
   dptr->deletedInfoIter->nxt = new infoList;
   dptr->deletedInfoIter      = dptr->deletedInfoIter->nxt;
   dptr->deletedInfoIter->inf = currInfo;
   dptr->deletedInfoIter->nxt = NULL;
};
void* frontpareto1::getInfoCurrElem(){
    return (*(dptr->mainIter))->info;
}
void frontpareto1::putInfoOnLastAdded(void* newinfo)
{
    dptr->lastAdded->info = newinfo;
}
#endif /*ATTACH_INFO_TO_PAIRS*/
