#ifndef __FILESPLIT_H_INCLUDED__
#define __FILESPLIT_H_INCLUDED__


#include <list>
#include <vector>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

class FileSplits 
{
	int splitSize; // in MB
	list<string> *contents;
	
	string buffer;
			
	public:
	FileSplits()
	{
	   contents = new list<string>();
	}
	
	~FileSplits() 
	{
	  contents->clear();
	  delete contents;
	}
	
	void setContents(list<string> *contents)
	{
	   this->contents = contents;
	}
	
	void clear()
	{
	  contents->clear();
	}
	
	list<string>* getContents() const;
	int numLines();
		
	void AddGeom(string geom); 
	void Write(char *buffer);
	void WriteV2(char *buffer);
	void WriteBuffer(char *buffer);
	void ChunkBuffer(int block_size);
	void PrintK(int k);
};

#endif 
