#include "driver.h"
#include "voro++.hh"
#include "math.h"

#define MAXLINE 512
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/*------------------------------------------------------------------------------
 * Method to prepare for the FEFF9 input for XANES or EXAFS calculations.
 *----------------------------------------------------------------------------*/
void Driver::FEFF_main()
{
  printf("\n"); for (int i=0; i<9; i++) printf("===="); printf("  FEFF  ");
  for (int i=0; i<9; i++) printf("===="); printf("\n");
  // selection of atoms for each frame
  char str[MAXLINE], select[MAXLINE], workdir[MAXLINE];
  one = all[0];

  while (1){
    printf("Please input your selection command, type h for help [all]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(select, str);
      char *ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0) one->SelHelp();
      else {
        one->selection(select);
        one->SelInfo();
      }
    } else break;
  }

  for (int i=0; i<20; i++) printf("===="); printf("\n");
return;
}
