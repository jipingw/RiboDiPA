#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector psiteCal(CharacterVector cigar, IntegerVector start, IntegerVector psitemap){
  
  int n=start.size();
  NumericVector center(n);
  
  /* prase cigar durations */
  for (int i=0; i < n; i++){
    /*Rcout << i<<std::endl;*/
    char temp[strlen(cigar[i])+1];
    strcpy(temp,cigar[i]);
    /*Rcout << "temp="<<temp<<std::endl;*/
    char * pch;
    pch= strtok(temp,"MIDNSHP=X");
    int j=0;
    
    while (pch != NULL)
    {
      /*Rcout << pch<<std::endl;*/
      pch = strtok (NULL, "MIDNSHP=X");
      j=j+1;
    }
    /*Rcout << "j="<<j<<std::endl;*/
    int n;
    n=j;
    if(n>3){
      Rcout << "  Warning: reads #" <<i+1 << " aligned to more than three consecutive regions! \n" << std::flush;
      R_FlushConsole();
      Rcout << "           The cigar string is:" << cigar[i] << "\n"<< std::flush;
      R_FlushConsole();
      Rcout << "           This may indicate this read is not a Ribo-seq read or alignment is invalid. \n"<< std::flush;
      R_FlushConsole();
    }
    
    NumericVector duration(n);
    j=0;
    strcpy(temp,cigar[i]);
    pch= strtok(temp,"MIDNSHP=X");
    
    while (pch != NULL)
    {
      duration[j]=std::atoi(pch);
      /*Rcout << pch<<std::endl;*/
      pch = strtok (NULL, "MIDNSHP=X");
      j=j+1;
    }
    duration=duration[Rcpp::Range(0,j-1)];
    center[i]= duration[0];
    /* prase letters*/
    strcpy(temp,cigar[i]);
    pch = strtok(temp, "0123456789");
    
    CharacterVector letter(n);
    j=0;
    while (pch != NULL)
    {
      letter[j]=pch;
      /*Rcout << pch<<std::endl;*/
      pch = strtok(NULL, "0123456789");
      j=j+1;
    }
    letter=letter[Rcpp::Range(0,j-1)];
    /*Rcout <<letter<<std::endl;*/  
    /*calculate read center*/
    int readlength=0;
    for (j=0;j<letter.size();j++){
      if(letter[j]=="M"||letter[j]=="I"||letter[j]=="S"||letter[j]=="="||letter[j]=="X"){
        readlength=readlength+duration[j];
      }
    }
    
    /*Rcout<<"read length "<<readlength<<std::endl;*/
    int fragmentlength=0;
    int position=start[i];
    
    if(duration.size()==1)
      center[i]=start[i]+psitemap[i];
    else{
      j=0;
      while(j<duration.size()){
        if(letter[j]=="M"||letter[j]=="X"||letter[j]=="="){
          fragmentlength=fragmentlength+duration[j];
          position=position+duration[j];
          if(fragmentlength>psitemap[i]){
            center[i]=position-(fragmentlength-psitemap[i]);
            j=duration.size();}
        }else if(letter[j]=="I")
          fragmentlength=fragmentlength+duration[j];
        else if(letter[j]=="S"){
          fragmentlength=fragmentlength+duration[j];
          if (j>0)
            position=position+duration[j];
        }else if(letter[j]=="N"||letter[j]=="D"){
          position=position+duration[j];
        }
        j=j+1;
      }
    }
    /*Rcout<<"fragmentlength "<<fragmentlength<<std::endl;*/
  }
  /*Rcout<<"center"<<center<<std::endl;*/
  return center;
}
