#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "./vcctl.h"

int Xsize,Ysize,Zsize;
int Xsizeorig,Ysizeorig,Zsizeorig;
float Version,Res;

int ***Mic;

/* Function declarations */
FILE *filehandler(char *prog, char *filename, char *tocheck);
int ***ibox(size_t xsize, size_t ysize, size_t zsize);
void free_ibox(int ***fc, size_t xsize, size_t ysize);
void bailout(char *name, char *msg);
int read_imgheader(FILE *fpin, float *ver, int *xsize,
	int *ysize, int *zsize, float *res);

int main(int argc, char *argv[])
{
    register int i,j,k;
    int ival,val;
    char buff[128],filename[128];
    FILE *fpin,*fpout;

    if (argc != 2) {
       printf("\n\nUSAGE:  vcctl2thames [image file name]\n\n");
       return(0);
    }

    fpin = filehandler("vcctl2thames",argv[1],"READ");

    /* Read in the vcctl microstructure image header */
    if (read_imgheader(fpin,&Version,&Xsize,&Ysize,&Zsize,&Res)) {
        fclose(fpin);
        bailout("vcctl2thames","Error reading image header");
        return(1);
    }

    Mic = ibox(Xsize,Ysize,Zsize);

    sprintf(filename,"%s.init",argv[1]);

    if ((fpout = fopen(filename,"w")) == NULL) {
        printf("\nERROR:  Could not open file %s. Exiting.\n\n",filename);
        free_ibox(Mic,Xsize,Ysize);
        exit(1);
    }

    fprintf(fpout,"Version: 5.0\n");
    fprintf(fpout,"X_Size: %d\n",Xsize);
    fprintf(fpout,"Y_Size: %d\n",Ysize);
    fprintf(fpout,"Z_Size: %d\n",Zsize);
    fprintf(fpout,"Image_Resolution: 1.0\n");

    /* Read in the vcctl microstructure image */

    for (k = 0; k < Zsize; k++) {
        for (j = 0; j < Ysize; j++) {
            for (i = 0; i < Xsize; i++) {
                fscanf(fpin,"%s",buff);
                ival = atoi(buff);
                switch (ival) {
                    case POROSITY:
                        val = 1;
                        break;
                    case C3S:
                        val = 2;
                        break;
                    case C2S:
                        val = 3;
                        break;
                    case C3A:
                        val = 4;
                        break;
                    case C4AF:
                        val = 5;
                        break;
                    case CAS2:
                        val = 6;
                        break;
                    case K2SO4:
                        val = 7;
                        break;
                    case NA2SO4:
                        val = 8;
                        break;
                    case GYPSUM:
                        val = 9;
                        break;
                    case HEMIHYD:
                        val = 10;
                        break;
                    case ANHYDRITE:
                        val = 10;
                        break;
                    case CACO3:
                        val = 11;
                        break;
                    case CH:
                        val = 12;
                        break;
                    case CSH:
                        val = 13;
                        break;
                    case AFMC:
                        val = 14;
                        break;
                    case AFM:
                        val = 15;
                        break;
                    case ETTR:
                        val = 16;
                        break;
                    case BRUCITE:
                        val = 17;
                        break;
                    case FREELIME:
                        val = 18;
                        break;
                    default:
                        val = 0;
                        break;
                }
                Mic[i][j][k] = val;
            }
        }
    }

    fclose(fpin);

    /* Now write the rest of the thames image file */

    for (k = 0; k < Zsize; k++) {
        for (j = 0; j < Ysize; j++) {
            for (i = 0; i < Xsize; i++) {
                fprintf(fpout,"%d\n",Mic[i][j][k]);
            }
        }
    }

    fclose(fpout);
    free_ibox(Mic,Xsize,Ysize);
    exit(0);
}

/* Function definitions */

FILE *filehandler(char *prog, char *filename, char *tocheck)
{
	FILE *fptr;

	fptr = NULL;

	if (!strcmp(tocheck,"NOCLOBBER")) {
		if ((fptr = fopen(filename,"r")) != NULL) {
			printf("\nERROR in %s:",prog);
			printf("\n\tFile %s already exists.",filename);
			printf("\n\tPlease verify file name. Program is ");
			printf("exiting now.\n\n");
			fflush(stdout);
			fclose(fptr);
			fptr = NULL;
		} else if ((fptr = fopen(filename,"w")) == NULL) {
			printf("\nERROR in %s:",prog);
			printf("\n\tCould not create file %s",filename);
			printf("\n\tPlease verify write permissions. Program is ");
			printf("exiting now.\n\n");
			fflush(stdout);
		}
	} else if (!strcmp(tocheck,"READ")) {

		if ((fptr = fopen(filename,"r")) == NULL) {
			printf("\nERROR in %s:",prog);
			printf("\n\tFile %s could not be opened for ",filename);
			printf("reading.");
			printf("\n\tPlease verify file name. Program is ");
			printf("exiting now.\n\n");
			fflush(stdout);
		}

	} else if (!strcmp(tocheck,"READ_NOFAIL")) {

		fptr = fopen(filename,"r");

	} else if (!strcmp(tocheck,"WRITE")) {

		if ((fptr = fopen(filename,"w")) == NULL) {
			printf("\nERROR in %s:",prog);
			printf("\n\tFile %s could not be created.",filename);
			printf("\n\tPlease verify file name. Program is ");
			printf("exiting now.\n\n");
			fflush(stdout);
		}

	} else if (!strcmp(tocheck,"APPEND")) {

		if ((fptr = fopen(filename,"a")) == NULL) {
			printf("\nERROR in %s:",prog);
			printf("\n\tFile %s could not be opened for ",filename);
			printf("appending.");
			printf("\n\tPlease verify file name. Program ");
			printf("is exiting now.\n\n");
			fflush(stdout);
		}

	}

	return(fptr);

}

/***
*	ibox
*
*	Routine to allocate memory for an 3D array of ints
*	All array indices are assumed to start with zero.
*
*	Arguments:	int number of elements in each dimension
*	Returns:	Pointer to memory location of first element
*
*	Calls:		no other routines
*	Called by:	main routine
*
***/
int ***ibox(size_t xsize, size_t ysize, size_t zsize)
{
	size_t i,j;
	int ***fc;

	fc = (int ***)malloc(xsize * sizeof(*fc));
	if (!fc) {
		printf("\n\nCould not allocate space for column of ibox.");
		return(NULL);
	}

    for (i = 0; i < xsize; ++i) {
        fc[i] = NULL;
    }

    for (i = 0; i < xsize; ++i) {
	    fc[i] = (int **)malloc(ysize * sizeof(*fc[i]));
        if (!fc[i]) {
		    printf("\n\nCould not allocate space for row of ibox.");
            free_ibox(fc,xsize,ysize);
		    return(NULL);
        }
    }

    for (i = 0; i < xsize; ++i) {
        for (j = 0; j < ysize; ++j) {
            fc[i][j] = NULL;
        }
    }

    for (i = 0; i < xsize; ++i) {
        for (j = 0; j < ysize; ++j) {
	        fc[i][j] = (int *)malloc(zsize * sizeof(*fc[i][j]));
            if (!fc[i][j]) {
		        printf("\n\nCould not allocate space for depth of ibox.");
                free_ibox(fc,xsize,ysize);
		        return(NULL);
            }
        }
    }

	return(fc);
}

/***
*	free_ibox
*
*	Routine to free the allocated memory for a 3D array of ints
*	All array indices are assumed to start with zero.
*
*	Arguments:	Pointer to memory location of first element
*	            x dimension of the box
*	            y dimension of the box
*
* 	Returns:	Nothing
*
*	Calls:		no other routines
*	Called by:	main routine
*
***/
void free_ibox(int ***fc, size_t xsize, size_t ysize)
{
    size_t i,j;
    if (fc != NULL) {
        for (i = 0; i < xsize; ++i) {
            if (fc[i] != NULL) {
                for (j = 0; j < ysize; ++j) {
                    if (fc[i][j] != NULL) free(fc[i][j]);
                }
                if (fc[i] != NULL) free(fc[i]);
            }
        }
    }
    fc = NULL;

	return;
}

/******************************************************************************
*	Function bailout prints error message to stdout
*
* 	Arguments:	char string of program name
* 				char error message
*
* 	Returns:	nothing
*
*	Programmer:	Jeffrey W. Bullard
*				NIST
*				100 Bureau Drive, Stop 8615
*				Gaithersburg, Maryland  20899-8615
*				USA
*
*				Phone:	301.975.5725
*				Fax:	301.990.6891
*				bullard@nist.gov
*				
*	16 March 2004
******************************************************************************/
void bailout(char *name, char *msg)
{
	printf("\nERROR in %s:",name);
	printf("\n\t%s",msg);
	printf("\n\tExiting now.\n");
	fflush(stdout);
	return;
}

/******************************************************************************
*    Function read_imgheader reads all the header information of
*    an image file (already assumed to be open), and then returns control
*    to calling function	for further reading of the microstructure itself
*
*     Arguments:    file pointer
*                   pointer to float version
*                   pointer to int xsize
*                   pointer to int ysize
*                   pointer to int zsize
*                   pointer to float resolution
*
*    Returns:    int status flag (0 if okay, 1 if otherwise)
*
*    Programmer:    Jeffrey W. Bullard
*                   NIST
*                   100 Bureau Drive, Stop 8615
*                   Gaithersburg, Maryland  20899-8615
*                   USA
*
*                   Phone:    301.975.5725
*                     Fax:    301.990.6891
*                             bullard@nist.gov
*				
*    16 March 2004
******************************************************************************/
int read_imgheader(FILE *fpin, float *ver, int *xsize,
	int *ysize, int *zsize, float *res)
{
    int status = 0;
    char buff[MAXSTRING],buff1[MAXSTRING];

    if (!fpin) {
        status = 1;
        return(status);
    }
		
    fscanf(fpin,"%s",buff);
    if (!strcmp(buff,VERSIONSTRING)) {
        fscanf(fpin,"%s",buff);
        *ver = atof(buff);
        fscanf(fpin,"%s",buff);
        if (!strcmp(buff,XSIZESTRING)) {
            fscanf(fpin,"%s",buff);
            *xsize = atoi(buff);
            fscanf(fpin,"%s %s",buff,buff1);
            *ysize = atoi(buff1);
            fscanf(fpin,"%s %s",buff,buff1);
            *zsize = atoi(buff1);
            fscanf(fpin,"%s %s",buff,buff1);
            *res = atof(buff1);
        } else if (!strcmp(buff,IMGSIZESTRING)){
            fscanf(fpin,"%s",buff);
            *xsize = atoi(buff);
            *ysize = *xsize;
            *zsize = *xsize;
            *res = 1.0;
       }

    } else {

        /***
        *    This image file was generated prior to
        *    Version 3.0.  Allow backward compatibility
        *    by defaulting system size to 100 and
        *    system resolution to 1.0
        ***/

        *ver = 2.0;
        *res = DEFAULTRESOLUTION;
        *xsize = DEFAULTSYSTEMSIZE;
        *ysize = DEFAULTSYSTEMSIZE;
        *zsize = DEFAULTSYSTEMSIZE;
	
        rewind(fpin);
    }
		
    return(status);
}
