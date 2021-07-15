#include "driver.h"
#include "timer.h"

/*------------------------------------------------------------------------------
 * Method to count the number of atoms within selection as a function of
 * MD steps (for selected frames)
 *----------------------------------------------------------------------------*/
void Driver::count_selected()
{
  char str[MAXLINE], *ptr, selcmd[MAXLINE];
  one = all[istr];
  while (1){
    printf("\nPlease input the atom selection command, `h` for help [all]: ");
    if (count_words(fgets(str,MAXLINE,stdin)) > 0){
      strcpy(selcmd, str);
      ptr = strtok(str," \n\t\r\f");
      if (strcmp(ptr,"h") == 0){ one->SelHelp(); continue; }
    } else strcpy(selcmd,"all\n");
  
    // check the selection command on the first frame
    one->selection(selcmd); one->SelInfo();
    if (one->nsel < 1){
      printf("It seems that no atom is selected, are you sure about this? (y/n)[y]: ");
      if (count_words(fgets(str,MAXLINE,stdin)) > 0){
        char *ptr = strtok(str," \n\t\r\f");
        if (strcmp(ptr,"y")!= 0 && strcmp(ptr,"Y")!=0) continue;
      }
    }
    break;
  }

  // output the result
  printf("\nPlease input the file name to output the result [count.dat]: ");
  fgets(str,MAXLINE, stdin);
  ptr = strtok(str, " \n\t\r\f");
  if (ptr == NULL){
    strcpy(str, "count.dat");
    ptr = strtok(str, " \n\t\r\f");
  }
  FILE *fp = fopen(ptr,"w");
  fprintf(fp,"# Number of atoms within selection wrt MD time: %s", selcmd);
  fprintf(fp,"# Istep  Number-atoms-in-selection\n");

  int nused = 0;
  // loop over all frames
  for (int img = istr; img <= iend; img += inc){
    one = all[img];

    // select atoms as source
    one->selection(selcmd);

    fprintf(fp, "%d %d\n", one->tstep, one->nsel);
    ++nused;
  }
  fclose(fp);

  printf("\n%d images were analyzed and the result is written to %s\n", nused, ptr);
  for (int i = 0; i < 20; ++i) printf("===="); printf("\n");

return;
}
/*------------------------------------------------------------------------------*/
