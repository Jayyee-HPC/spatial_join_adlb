#include "FileSplits.h"

list<string>* FileSplits :: getContents() const
{
  return contents;
}

void FileSplits :: AddGeom(string geom) 
{
  contents->push_back(geom);
}


int FileSplits :: numLines()
{
  if(contents!=nullptr)
  {
    return contents->size();
  }
  else
  {
    return 0;
  }
}

void FileSplits :: PrintK(int k)
{
 int counter = 0;
 list<string>::const_iterator i;

 for(i = contents->begin(); i!= contents->end() && counter < k; i++) {
    cout<<*i<<endl;
    counter++;
  }
}

void FileSplits :: WriteV2(char *data)
{
   std::stringstream ss(data);
   std::string to;
 size_t total = 0;
  if (data != NULL)
  {
    while(std::getline(ss,to,'\n')){
      //cout << to <<endl;
      contents->push_back(to);
      total++;
    }
  }

  cout<<" FileSplits total "<<total<<endl;

}


void FileSplits :: Write(char *data){
    // to avoid modifying original string
    // first duplicate the original string and return a char pointer then free the memory
    char *dup = strdup(data);
    
    char *token = std::strtok(dup, "\n");
    while(token != NULL){
        contents->push_back(string(token));
        // the call is treated as a subsequent calls to strtok:
        // the function continues from where it left in previous invocation
        token = std::strtok(NULL, "\n");
    }
    free(dup);
}

void FileSplits :: WriteBuffer(char *data){
    // to avoid modifying original string
    // first duplicate the original string and return a char pointer then free the memory
    buffer = strdup(data);
}

void FileSplits :: ChunkBuffer(int block_size) {
	
	vector<int> pos;
	pos.push_back(-1);
	
	int i, j=1;
	for(i=0; buffer[i] != '\0'; i++){
		if(buffer[i] == '\n')
			j++;
		
		if (j%block_size == 0) {
			pos.push_back(i);
			j=1;
		}
	}
	if (*(pos.end()) != i-1)
		//cout <<"FileSplits end pointer: " << *(pos.end()) << " i-1:"<< i - 1<<endl;
		pos.push_back(i-1);
	
	//long charCount = 0;
	for(int k=0; k< ((int)pos.size())-1; k++) {
		string sub_str = buffer.substr(pos[k]+1,(pos[k+1] - pos[k]));
		
		//charCount += sub_str.size();

		contents->push_back(sub_str);
	}
	//printf("file part %d, file length %lld \n", i, charCount);
	
}

