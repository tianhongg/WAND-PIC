//----------------------------------------------------------------------------------||
//----------  This Piece Of Code Is Borrowed From The PiC Simulation Code:  --------||
//----------  Virtual Laser Plasma Lab (VLPL). The Author Of Wand-PIC Would --------||
//----------  Like To Thank Whoever Wrote This Code.                        --------||
//----------------------------------------------------------------------------------||



#ifndef H_NLIST
#define H_NLIST
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#define BUF_SIZE 16384


class NList_Entry;

class NList {
private:
   char* list_name;
   NList_Entry* head;
   NList_Entry* tail;
   char buf[BUF_SIZE];
   int bsize;
   static char start_tag; // Character to mark start of NList in 
   
   static char close_tag; // Character to mark end of NList in file
   static char rem_tag;   // Character to mark comment in NList file
   void error(int ierr);
   char* err_msg[10];
public:

   NList_Entry* AddEntry(char* n, int* ptr, int v = 0);      //  a function under the class NList_Entry
   NList_Entry* AddEntry(char* n, long* ptr, long v = 0);
   NList_Entry* AddEntry(char* n, float* ptr, float v = 0.);
   NList_Entry* AddEntry(char* n, double* ptr, double v = 0.);
   NList_Entry* AddEntry(char* n, char* ptr, char* v =(char *)"%s");

   int read(FILE* f); // reads the file f containing NList, returns 0 if OK;
   int sread(char* fc); // reads from string fc containing NList, returns 0 if OK;

   int write(FILE* f); // writes to file f, returns 0 if OK;
   int swrite(char* fc); // writes to string fc, returns 0 if OK;

   NList(const char* name);
   NList();
   virtual ~NList(){}
   
   friend class NList_Entry;
};

class NList_Entry {
private:
   NList* list; // The Namelist the entry belongs to
   char* name;  // name of the variable
   void* ptr;   // pointer to the variable itself
   char* fmt;   // format the variable must be converted, like "%g"
   char t;
   NList_Entry* next; // the next entry;
   int size;

public:
   void setsize(int a_size) { size = a_size;};

private:
   NList_Entry(NList* l, char* n, int* p, int v);
   NList_Entry(NList* l, char* n, long* p, long v);
   NList_Entry(NList* l, char* n, float* p, float v);
   NList_Entry(NList* l, char* n, double* p, double v);
   NList_Entry(NList* l, char* n, char* p, char* v);

   friend class NList;
};

#endif
