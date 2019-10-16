//----------------------------------------------------------------------------------||
//----------  This Piece Of Code Is Borrowed From The PiC Simulation Code:  --------||
//----------  Virtual Laser Plasma Lab (VLPL). The Author Of Wand-PIC Would --------||
//----------  Like To Thank Whoever Wrote This Code.                        --------||
//----------------------------------------------------------------------------------||


#include "namelist.h"


char NList::start_tag = '&';
char NList::close_tag = '@';
char NList::rem_tag = '#';

// --------------- Int   Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, int* ptr, int v)
{
   NList* l = this;   // point to this function
   NList_Entry* e = new NList_Entry(l, n, ptr, v); // point to this obj
   e->next=head; head = e;
   return e;
}

// --------------- Long   Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, long* ptr, long v)
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- Float Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, float* ptr, float v) 
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- Double Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, double* ptr, double v) 
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- String Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, char* ptr, char* v) 
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- NList::NList ----------------------------------------

NList::NList() { // This constructor has been made for Controls::C_f_SaveTime.
//   NList::start_tag = '&';
//   NList::close_tag = '/';
//   NList::rem_tag = '#';

   tail=head=NULL;

   err_msg[0] = (char*)"FILE_NOT_OPEN";
   err_msg[1] = (char*)"NLIST_NOT_FOUND";
   err_msg[2] = (char*)"WRONG_MEMBER";
   err_msg[3] = (char*)"CONVERSION_ERROR";
   err_msg[4] = (char*)"PREMATURE_END_OF_FILE";
}


// --------------- NList::NList ----------------------------------------

NList::NList(const char* name) {
   list_name = new char[strlen(name)+2];
   sprintf(list_name,"%s",name);
   tail=head=NULL;
   bsize=BUF_SIZE;
   err_msg[0] = (char*)"FILE_NOT_OPEN";
   err_msg[1] = (char*)"NLIST_NOT_FOUND";
   err_msg[2] = (char*)"WRONG_MEMBER";
   err_msg[3] = (char*)"CONVERSION_ERROR";
   err_msg[4] = (char*)"PREMATURE_END_OF_FILE";
}

// --------------- Read NList ----------------------------------------

int NList::read(FILE* f) 
{
char *stmp = NULL; 
   int found = -1;
   //  cout << "Looking for Namelist:" << list_name <<endl;

   if (f == NULL) error(1);
   while (fgets(buf,bsize,f)) {


      if (buf[0] == rem_tag) continue;
      if ((stmp=strchr(buf,start_tag))) {
         stmp = strtok(stmp+1," \t\n");
         if (strcmp(stmp,list_name)) continue;
         break; 
      };
   };

   if (stmp==NULL) error(2);
   while ((stmp=fgets(buf,bsize,f))) {
      if (buf[0] == rem_tag) continue;
      if (strchr(buf,close_tag)) {
         //      cout << "Namelist:" << list_name << " read" <<endl;
         return 0;
      }
      if ((stmp = strtok(stmp," =\t\n")) == NULL) continue;
      if (strcmp(stmp,"\n")==0) continue;
      NList_Entry* e = head;
      while (e) {
         if (strcmp(stmp,e->name)) { e=e->next; continue; }
         if ((stmp = strtok(NULL," =\t\n")) == NULL) continue;
         switch (e->t) {
      case 'i':  if (sscanf(stmp,"%d",(int*)e->ptr)==0) error(4); break;
      case 'l':  if (sscanf(stmp,"%ld",(long*)e->ptr)==0) error(4); break;
      case 'f':  if (sscanf(stmp,"%g",(float*)e->ptr)==0) error(4); break;
      case 'd':  
         { float dtmp; double* dptr;
         if (sscanf(stmp,"%g",&dtmp)==0) error(4);
         dptr=(double*)e->ptr; *dptr = dtmp; 
         break;
         }
      case 's':  if (sscanf(stmp,e->fmt,(char*)e->ptr)==0) error(4); 
         break;
      default: 
         ;}
         break;
      }
      if (stmp==NULL) continue;
      if (e==NULL) continue; //<SergK>{ cout << stmp << "\n"; error(3); }
   }

   error(5);
   return -5;



}


// --------------- NList Error Handler ------------------------------------

void NList::error(int ierr) 
{
   std::cout << "NList "<< list_name<<" error: "<<err_msg[ierr-1] <<'\n';
   exit(ierr);
}

// --------------- Write NList ----------------------------------------


// --------------- Read NList ----------------------------------------


// --------------- NList Constructors ------------------------------------

NList_Entry::NList_Entry(NList* l, char* n, char* p, char* v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   fmt = v;
   t = 's';
   int m_size = 0;
   if (sscanf(v+1,"%d",&size)) {       
      m_size = size;
   } else {
      size = 1;
   };
}

NList_Entry::NList_Entry(NList* l, char* n, int* p, int v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   *p = v;
   fmt = (char*)"%d";
   t = 'i';
   size = sizeof(int);
}

NList_Entry::NList_Entry(NList* l, char* n, long* p, long v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   *p = v;
   fmt = (char*)"%d";
   t = 'l';
   size = sizeof(long);
}

NList_Entry::NList_Entry(NList* l, char* n, float* p, float v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   *p = (float)v;
   fmt = (char*)"%g";
   t = 'f';
   size = sizeof(float);
}

NList_Entry::NList_Entry(NList* l, char* n, double* p, double v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   *p = (double)v;
   fmt = (char*)"%g";
   t = 'd';
   size = sizeof(double);
}
