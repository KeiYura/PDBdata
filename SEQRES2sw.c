/******************************************************************************/
/*                                                                            */
/*     SEQRES2sw.c                                                            */
/*               Convert PDB SEQRES line  to SwissProt format                 */
/*                                                                            */
/*                                 Kei Yura                                   */
/*                                             Feb. 19th, 2002                */
/*                                             Jul. 19th, 2003                */
/*                                             Jul. 20th, 2003                */
/*                                             Aug. 19th, 2003                */
/*                                             May  24th, 2006                */
/*                                                                            */
/******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#define MAXLINE 100
#define MAXRES  3000
#define MAXCH   62


int LIMIT;
char *One;
char **Three;
FILE *err;
int   I;
int   CHAIN;
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void init()
{
	int  total, i;
	char file[MAXLINE];
	char lines[MAXLINE];
	FILE *fp;
	
	file[0] = '\0';
	strcat(file,getenv("HOME"));
	strcat(file,"/parm/residue.dat");

	if ((fp = fopen(file,"r")) == NULL)
	{
	     fprintf(stderr,"Cannot open %s\n",file);
	     exit (-1);
	}
	
	if ((err = fopen("SEQRES2sw.log","a")) == NULL)
	{
	     fprintf(stderr,"Cannot open ATOMtosw.log\n");
	     exit (-1);
	}

	total = 0;
	while(fgets(lines,MAXLINE,fp) != NULL)
	     total++;
	
	if ((One = (char *)(malloc(total*sizeof(char)))) == NULL)
	{
	     fprintf(stderr,"Cannot allocate memory for *One\n");
	     exit (-1);
	}

	if ((Three = (char **)(malloc(total*sizeof(char *)))) == NULL)
	{
	     fprintf(stderr,"Cannot allocate memory for *Three\n");
	     exit (-1);
	}

	for (i = 0; i < total; i++)
	{
	     if ((Three[i] = (char *)(malloc(4*sizeof(char)))) == NULL)
	     {
		  fprintf(stderr,"Cannot allocate memory for Three[]\n");
		  exit (-1);
	     }
	}

	rewind(fp);

	total = 0;
	while(fgets(lines,MAXLINE,fp) != NULL)
	{
	     One[total] = lines[4];
	     strncpy(Three[total],lines,3);
	     Three[total][3] = '\0';
	     total++;
	}

	One[total] = '\0';

	return;
}

/******************************************************************************/
/******************************************************************************/
char t_to_o(s)
char *s;
{
        int i, len;

        len = strlen(One);
        for (i = 0; i < len; i++)
             if(strncmp(s,Three[i],3) == 0)
                  return(One[i]);

        return('X');
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void conv(fp)
FILE *fp;
{
	char    lines[MAXLINE], id[MAXLINE], chain[MAXCH];
	char    dt[MAXLINE], de[MAXCH][MAXLINE];
        char    os[MAXCH][MAXLINE], oc[MAXCH][MAXLINE];
	char    *seq, ch[2], sq[MAXLINE], idp[MAXLINE], idp2[MAXLINE];
	char    dummy[MAXLINE], oldch[MAXCH];
	time_t  secs_now;
	int     molid[MAXCH];
	int     chnum = 0, chflag = 0, inside, cline = 0;
	int     i, j, row, nres, source = -1, compnd = -1, len, flag;
	void    makeseq();
	void    printseq();

	I = 0;
	CHAIN = 0;
	id[0] = '\0';
	while(fgets(lines,MAXLINE,fp) != 0)
	{
	   if (strncmp(lines,"HEADER ",7) == 0)
	   {                 /* HEADER is converted to ID */

	      len = strlen(lines);
	      for (i = 0; 62+i < len && lines[62+i] != ' ' 
                                  && lines[62+i] != '\n'; i++)
	     id[i] = lines[62+i];
	      id[i] = '\0';
	      dummy[0] = '\0';

	      if (strlen(id) == 0)
	      {
	           strcat(id,"XXXX");
                   fprintf(err,"Warning: NO ID exists\n");
	      }

	      time(&secs_now);
	      dt[0] = '\0';
	      sprintf(dt,"DT   converted from PDB (SEQRES  ) %s :%s",
	  		                         id, ctime(&secs_now));
	   }
           else if (strncmp(lines,"EXPDTA    THEORETICAL MODEL",27) == 0)
           {
              fprintf(err,"Warning: %s is a theoretical model\n",id);
           }
	   else if (strncmp(lines,"COMPND ",    7) == 0)
           {                           /* COMPND is converted to DE */
	      if (strncmp(&lines[10],"MOL_ID:",7) == 0 ||
	          strncmp(&lines[11],"MOL_ID:",7) == 0)
	      {
		   compnd++;
	           molid[compnd] = atoi(&lines[18]);
	           de[compnd][0] = '\0';
	      }
	      else if (strncmp(&lines[11],"MOLECULE: ",10) == 0)
	      {
	           len = strlen(lines);
	           for (i = 0 ; 21+i < len && i < 70 && lines[21+i] != ';' &&
                                          lines[21+i] != '\n'; i++)
		        dummy[i] = tolower(lines[21+i]);
	           dummy[i] = '\0';

	           if (dummy[0] == ' ' && dummy[1] == ' ')
	           {
	                if (fgets(lines,MAXLINE,fp) == NULL)
	                {
			     fprintf(stderr,"File format error %s\n",id);
			     exit (-1);
	                }
	 
	                for (i = 0 ; 11+i < len && i < 70 && lines[11+i] != ';' 
                                            && lines[11+i] != '\n'; i++)
		             dummy[i] = tolower(lines[11+i]);
	                dummy[i] = '\0';
	           }

	           sprintf(de[compnd],"DE   %s\n",dummy);
	      }
	      else if (strncmp(&lines[11],"CHAIN:",6) == 0)
	      {
	            if (strncmp(&lines[18],"NULL",4) == 0)
	                 chain[compnd] = ' ';
	            else
	            {
	                 len = strlen(lines);
	                 for (i = 18, j = 0; i < len && 
                           i < strlen(lines)-1 && lines[i] != ' '; i+=3,j++)
	                      chain[compnd+j] = lines[i];

	                 if (j > 1)
	                      for (i = 1; i < j; i++)
			      {
			           strcpy(de[compnd+i],de[compnd]);
			           molid[compnd+i] = molid[compnd];
	                      }
		         compnd += (j-1);
		    }
	         
		    if (compnd >= MAXCH)
		    {
		         fprintf(stderr,"increase MAXCH (%s).\n",id);
			 exit (-1);
	            }
	      }
	      else if (strncmp(&lines[11],"EC:",3) == 0)
	           ;
	      else if (strncmp(&lines[11],"BIOLOGICAL_UNIT:",16) == 0)
	           ;
	      else if (strncmp(&lines[11],"FRAGMENT:",9) == 0)
	           ;
	      else if (strncmp(&lines[11],"SYNONYM:",8) == 0)
	           ;
	      else if (strncmp(&lines[11],"ENGINEERED:",11) == 0)
	           ;
	      else
	      {
	           len = strlen(lines);
		   if (compnd == -1)
		   {
	                for (i = 0 ; i < 60 && 10+i < len-1 
                                            && lines[10+i] != '\n'; i++)
		             dummy[i] = tolower(lines[10+i]);
	                dummy[i] = '\0';
			compnd++;
	                sprintf(de[compnd],"DE   %s\n",dummy);
		   	molid[compnd] = 1;
		        chain[compnd] = ' ';
	           }
	           else
	           {
		        if (strlen(de[compnd]) == 0)
			{

                             for (i = 0 ; i < 50 && lines[11+i] != ';' && 
                                      11+i < len && lines[11+i] != '\n'; i++)
                                  dummy[i] = tolower(lines[11+i]);
                             dummy[i] = '\0';
                             sprintf(de[compnd],"DE   %s\n",dummy);
	                }
	           }
	      }
	   }
	   else if (strncmp(lines,"SOURCE ",    7) == 0)
	   {
	      if (strncmp(&lines[10],"MOL_ID:",7) == 0 ||
	          strncmp(&lines[11],"MOL_ID:",7) == 0)
	      {
	                source = atoi(&lines[18]);
	      }
	      else if (strncmp(&lines[11],"ORGANISM_SCIENTIFIC: ",21) == 0)
	      {                            /* SOURCE is converted to OS */
	           len = strlen(lines);
	           for (i = 32, j = 0; i < len-1 && lines[i] != ';'; i++)
	           {
		        if ((lines[i] >= 'A' && lines[i] <= 'Z') ||
		            (lines[i] >= 'a' && lines[i] <= 'z') ||
		            (lines[i] >= '0' && lines[i] <= '9') ||
		             lines[i] == ' ' || lines[i] <= '.')
	                {
	                     if (j == 0)
		                  dummy[j++] = toupper(lines[i]);
	                     else
		                  dummy[j++] = tolower(lines[i]);
	                }
	           }
		   dummy[j] = '\0';

		   for (i = 0; i < compnd+1; i++)
	           {
		        if (molid[i] == source)
	                {
	                     sprintf(os[i],"OS   %s\n",dummy);
	                     sprintf(oc[i],"OC   %s\n",dummy);
	                }
	           }
	     }
	     else if (strncmp(&lines[11],"ORGANISM_COMMON:",16) == 0)
             {                            /* SOURCE is converted to OS */
	           len = strlen(lines);
                   for (i = 28, j = 0;
                              i < len-1 && lines[i] != ';'; i++)
                        dummy[j++] = tolower(lines[i]);
                   dummy[j] = '\0';

                   for (i = 0; i < compnd+1; i++)
                        if (molid[i] == source)
	                {
                             sprintf(oc[i],"OC   %s\n",dummy);
	                }
             }
	     else if (strncmp(&lines[11],"SYNTHETIC: YES",14) == 0)
	     {
                   for (i = 0; i < compnd+1; i++)
                   {
                        if (molid[i] == source)
                        {
                             sprintf(os[i],"OS   chemically synthesized\n");
                             sprintf(oc[i],"OC   chemically synthesized\n");
                        }
                   }
             }
	     else if (strncmp(&lines[11],"EXPRESSION_SYSTEM:",18) == 0)
	          ;
	     else if (strncmp(&lines[11],"EXPRESSION_SYSTEM_COMMON:",25) == 0)
	          ;
	     else
	     {
	          if (source == -1)
	          {
	               len = strlen(lines);
                       for (i = 10, j = 0;
		           i < len-1 && !(lines[i-1] == ' ' &&
                                                         lines[i] == ' '); i++)
                       {
                            if ((lines[i] >= 'A' && lines[i] <= 'Z') ||
                                (lines[i] >= 'a' && lines[i] <= 'z') ||
                                (lines[i] >= '0' && lines[i] <= '9') ||
                                 lines[i] == ' ' || lines[i] == '.')
                            {
                                 if (j == 0)
                                      dummy[j++] = toupper(lines[i]);
                                 else
                                      dummy[j++] = tolower(lines[i]);
                            }
                       }
                       dummy[j] = '\0';
	               sprintf(os[0],"OS   %s\n",dummy);
	               sprintf(oc[0],"OC   %s\n",dummy);
	               source = -2;
	         }
	      }
	   }
	   else if (strncmp(lines,"SEQRES ",7) == 0)
	   {
	        inside = 0;
	        do 
	        {
	             if (inside == 0)
	             {
	                  makeseq(lines,ch,&nres,&seq,id,fp);
	                  inside = 1;
	             }

	             if (ch[0] == lines[11])
		     {
	                  for (i = 0; i < 13 && I < nres; i++)
	                       seq[I++] = t_to_o(&lines[19+(i*4)]);
	             }
		     else 
		     {
	                  if (nres >= LIMIT)
	                       printseq(seq,id,ch,nres,dt,de,os,oc,0);
	                  free((void *)seq);

	                  makeseq(lines,ch,&nres,&seq,id,fp);
                          inside = 1;
	                  for (i = 0; i < 13 && I < nres; i++)
	                       seq[I++] = t_to_o(&lines[19+(i*4)]);
	             }
	        } while (strncmp(fgets(lines,MAXLINE,fp),"SEQRES ",7) == 0);
	   }
	   else if (strncmp(lines,"ATOM ",5) == 0 ||
                    strncmp(lines,"END",3) == 0 ||
	            strncmp(lines,"ENDMDL",6) == 0)
	      break;
	}

	if (I != 0 && nres > LIMIT)
	     printseq(seq,id,ch,nres,dt,de,os,oc,0);
	
	return;
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int main(argc, argv)
int argc;
char *argv[];
{
	void init();
        void conv();

	if (argc > 1)
	     LIMIT = atoi(argv[1]);
	else
	     LIMIT = 0;

	init();
        conv(stdin);
	fclose(err);
        return 0;
}

/******************************************************************************/
/******************************************************************************/
void makeseq(lines,ch,nres,seq,id,fp)
char   lines[MAXLINE];
char   ch[2];
int   *nres;
char **seq;
char   id[MAXLINE];
FILE  *fp;
{
	I = 0;
        ch[0] = lines[11]; ch[1] = '\0';

        if (!((ch[0] >= 'a' && ch[0] <= 'z') ||
              (ch[0] >= 'A' && ch[0] <= 'Z') ||
              (ch[0] >= '0' && ch[0] <= '9') ||
               ch[0] == ' '))
        {
             fprintf(err, "Warning: %s has an irregular ",id);
             fprintf(err, "character in chain identifier: %c",ch[0]);
             fprintf(err, "\n          Ignored.\n");
             while(fgets(lines,MAXLINE,fp) != NULL)
	     {
                  if (strncmp(fgets(lines,MAXLINE,fp),"SEQRES ",7) == 0 &&
                     ((ch[0] >= 'a' && ch[0] <= 'z') ||
                      (ch[0] >= 'A' && ch[0] <= 'Z') ||
                      (ch[0] >= '0' && ch[0] <= '9') ||
                       ch[0] == ' ') )
                       break;
	     }
        }

        *nres = atoi(&lines[12]);
	if (*nres <= 0)
	     exit (-1);

        if ((*seq = (char *)( malloc(sizeof(char)*(*nres+1)) )) == NULL)
        {
             fprintf(stderr,"Malloc error (%s).\n",id);
             exit (1);
        }

	return;
}

/******************************************************************************/
/******************************************************************************/
void printseq(seq,id,ch,nres,dt,de,os,oc,cmpnd)
char *seq;
char  id[MAXLINE];
char  ch[2];
int   nres;
char  dt[MAXLINE];
char  de[MAXCH][MAXLINE];
char  os[MAXCH][MAXLINE];
char  oc[MAXCH][MAXLINE];
int   cmpnd;
{
	int  i, j;
	int  len;
	int  row;
	int  xnum;
	char idp[MAXLINE];
	char idp2[MAXLINE];
	char sq[MAXLINE];

 	seq[I] = '\0';

	xnum = 0;
	for (i = 0; i < nres; i++)
	     if (seq[i] == 'X') xnum++;
	if (xnum > nres/3)
	     return;

        sprintf(idp,"ID   %s%s_",id,ch);
        len = strlen(idp);

        for (i = len; i < 20; i++)
             strcat(idp," ");
        sprintf(idp2,"%sSTANDARD;      PRT;%6d AA.",idp,nres);
        printf("%s\n",idp2);
        printf("%s",dt);

        if (de[cmpnd][0] != '\0')
             printf("%s",de[cmpnd]);
        else
             printf("DE   unknown\n");

        if (strncmp(os[cmpnd],"OS ",3) == 0)
             printf("%s",os[cmpnd]);
        else
             printf("OS   unknown\n");

        if (strncmp(oc[cmpnd],"OC ",3) == 0)
             printf("%s",oc[cmpnd]);
        else
             printf("OC   unknown\n");

        printf("CC   This is generated from PDB SEQRES rows.       \n");
        printf("CC   The sequence is not always correct.         \n");
        sprintf(sq, "SQ   SEQUENCE%6d AA;      MW;          CN;\n",nres);

        printf("%s",sq);
        for (row = 0; row < (nres-1)/60+1; row++)
        {
             printf("     ");
             for (i = 0; i < 60 && row*60+i < nres; i++)
             {
                  printf("%c",seq[row*60+i]);
                  if ((int)((i+1)/10.0)*10 == i+1)
                       printf(" ");
             }
             printf("\n");
        }
        printf("//\n");

	return;
}
/* EOF */
