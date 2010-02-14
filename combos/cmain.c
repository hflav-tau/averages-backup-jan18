#include <stdlib.h>
#include <string.h>

void combos_();

int main(int narg, char **argv){
  
  char *file=" ";
  int  i,len=1;

  if (narg>1) {
    for (i=1;i<narg;i++){
      len = strlen(argv[i]);
      combos_(&len,argv[i],&len);
    }
  } else {
    combos_(&len,file,&len);
  }

  return 0;

}
